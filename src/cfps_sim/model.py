"""Reduced, uncalibrated resource model. Concentrations: mM; time: minutes."""
import antimony
import math
from types import SimpleNamespace
import numpy as np
import torch
from torchdiffeq import odeint
from .backend import cuda_device, metadata
from .construct import AMINO_ACIDS, demands

# Illustrative effective parameters, not extracted from or fitted to competition data.
# ATP use for charging is lumped into ADP equivalents; this is not detailed adenylate chemistry.
MODEL = """
model cfps()
  species ATP, ADP, GTP, GDP, CTP, UTP, AA, Fuel, mRNA, Protein, GFP, Machinery;
  nt_a=180; nt_g=180; nt_c=180; nt_u=180; aa_length=240;
  tx := k_tx*dna_template*Machinery*(ATP/(K_nt+ATP))*(GTP/(K_nt+GTP))*(CTP/(K_nt+CTP))*(UTP/(K_nt+UTP));
  tl := k_tl*mRNA*Machinery*(ATP/(K_nt+ATP))*(GTP/(K_nt+GTP))*(AA/(K_aa+AA));
  mRNA' = tx-k_mrna*mRNA;
  ATP' = regen*Fuel/(K_fuel+Fuel)*ADP/(K_nt+ADP)-nt_a*tx-2*aa_length*tl;
  ADP' = -regen*Fuel/(K_fuel+Fuel)*ADP/(K_nt+ADP)+2*aa_length*tl;
  GTP' = -nt_g*tx-2*aa_length*tl;
  GDP' = 2*aa_length*tl;
  CTP' = -nt_c*tx;
  UTP' = -nt_u*tx;
  AA' = -aa_length*tl;
  Fuel' = -regen*Fuel/(K_fuel+Fuel)*ADP/(K_nt+ADP);
  Protein' = tl-k_mat*Protein;
  GFP' = k_mat*Protein;
  Machinery' = -k_decay*Machinery;
  ATP=1.5; ADP=0.5; GTP=1.5; GDP=0; CTP=1.5; UTP=1.5;
  AA=40; Fuel=20; mRNA=0; Protein=0; GFP=0; Machinery=1;
  dna_template=0.000005; k_tx=0.01; k_tl=0.2; k_mrna=0.01; k_mat=0.03;
  k_decay=0.002; regen=0.05; K_nt=0.1; K_aa=1; K_fuel=1;
end
"""


PARAMETERS = {"dna_template", "k_tx", "k_tl", "k_mrna", "k_mat", "k_decay", "regen", "K_nt", "K_aa", "K_fuel"}
INITIAL_SPECIES = {"ATP", "ADP", "GTP", "GDP", "CTP", "UTP", "Fuel", "mRNA", "Protein", "GFP", "Machinery"}
PARAMETER_ORDER = ("dna_template", "k_tx", "k_tl", "k_mrna", "k_mat", "k_decay",
                   "regen", "K_nt", "K_aa", "K_fuel")


@torch.compile
def resource_rates(y, p, ntp_counts, aa_counts):
    """Fuse small elementwise flux calculations into GPU kernels."""
    pools = y.clamp_min(0)
    ntps = torch.stack((pools[0], pools[2], pools[4], pools[5]))
    ntp_saturation = ntps / (p[7] + ntps)
    tx = p[1] * p[0] * pools[10] * ntp_saturation.prod()
    aa = pools[11:]
    aa_saturation = aa / (p[8] + aa)
    limitation = torch.where(aa_counts > 0, aa_saturation, torch.inf).min()
    tl = p[2] * pools[7] * pools[10] * ntp_saturation[0] * ntp_saturation[1] * limitation
    regen = p[6] * pools[6] / (p[9]+pools[6]) * pools[1] / (p[7]+pools[1])
    energy = 2 * aa_counts.sum() * tl
    maturation = p[4] * pools[8]
    rates = torch.stack((regen-ntp_counts[0]*tx-energy, -regen+energy,
                         -ntp_counts[1]*tx-energy, energy,
                         -ntp_counts[2]*tx, -ntp_counts[3]*tx, -regen,
                         tx-p[3]*pools[7], tl-maturation, maturation,
                         -p[5]*pools[10]))
    return torch.cat((rates, -aa_counts*tl))


def antimony_source(construct=None, parameters=None, initial=None):
    model = MODEL
    derived = None
    if construct is not None:
        if construct["target"].lower() != "sfgfp":
            raise ValueError("Sequence simulation currently supports sfGFP only; enzyme folding/activity models are not implemented")
        derived = demands(construct)
        counts = derived["aa_per_protein"]
        model = model.replace("UTP, AA, Fuel", "UTP, Fuel")
        individual = [f"aa_{aa}" for aa in AMINO_ACIDS]
        model = model.replace("  nt_a=", "  species " + ", ".join(individual) + ";\n  nt_a=")
        model = model.replace("  AA' = -aa_length*tl;", "\n".join(f"  aa_{aa}' = -{counts[aa]}*tl;" for aa in AMINO_ACIDS))
        model = model.replace("AA=40;", " ".join(f"aa_{aa}=2;" for aa in AMINO_ACIDS) + "\n  AA := " + "+".join(individual) + ";")
        factors = [f"aa_{aa}/(K_aa+aa_{aa})" for aa in AMINO_ACIDS if counts[aa]]
        limiting = factors[0]
        for factor in factors[1:]:
            limiting = f"min({limiting},{factor})"
        model = model.replace("(AA/(K_aa+AA))", f"({limiting})")
        for base, count in derived["ntp_per_transcript"].items():
            model = model.replace(f"nt_{base.lower()}=180", f"nt_{base.lower()}={count}")
        model = model.replace("aa_length=240", f"aa_length={derived['protein_aa']}")
    allowed_initial = INITIAL_SPECIES | ({f"aa_{aa}" for aa in AMINO_ACIDS} if construct else {"AA"})
    for values, allowed, kind in [(parameters or {}, PARAMETERS, "parameter"), (initial or {}, allowed_initial, "initial species")]:
        for name, value in values.items():
            if name not in allowed:
                raise ValueError(f"Unknown {kind}: {name}; chemistry absent from this model cannot be inferred")
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
            if name.startswith("K_") and value == 0:
                raise ValueError("Saturation constants must be positive")
            model = model.replace("end\n", f"  {name}={value:.17g};\nend\n")
    return model


class Trajectory(np.ndarray):
    """Host export of a completed GPU trajectory, with named columns."""
    def __new__(cls, values, colnames):
        result = np.asarray(values).view(cls)
        result.colnames = list(colnames)
        return result

    def __array_finalize__(self, obj):
        self.colnames = getattr(obj, "colnames", None)


class GPUModel:
    """Sequence-aware ODEs on CUDA; Antimony is used only for SBML export."""
    def __init__(self, construct=None, parameters=None, initial=None):
        # Validate inputs before allocating device state.
        self.source = antimony_source(construct, parameters, initial)
        self.device = cuda_device()
        self.integrator = SimpleNamespace(absolute_tolerance=1e-14, relative_tolerance=1e-10)
        self.parameters = dict(dna_template=0.000005, k_tx=0.01, k_tl=0.2,
                               k_mrna=0.01, k_mat=0.03, k_decay=0.002,
                               regen=0.05, K_nt=0.1, K_aa=1., K_fuel=1.)
        self.parameters.update(parameters or {})
        derived = demands(construct) if construct else None
        self.counts = derived["ntp_per_transcript"] if derived else dict(A=180, G=180, C=180, U=180)
        self.length = derived["protein_aa"] if derived else 240
        self.aa_counts = ([derived["aa_per_protein"][aa] for aa in AMINO_ACIDS]
                          if derived else [240])
        self.species = ["ATP", "ADP", "GTP", "GDP", "CTP", "UTP", "Fuel",
                        "mRNA", "Protein", "GFP", "Machinery"]
        self.species += [f"aa_{aa}" for aa in AMINO_ACIDS] if construct else ["AA"]
        pools = [1.5, 0.5, 1.5, 0., 1.5, 1.5, 20., 0., 0., 0., 1.]
        pools += [2.] * 20 if construct else [40.]
        self.initial = dict(zip(self.species, pools))
        self.initial.update(initial or {})
        self.original_parameters = self.parameters.copy()
        self.original_initial = self.initial.copy()
        self.state = self.initial.copy()
        self.construct = construct
        self.stoichiometry = torch.tensor(self.aa_counts, dtype=torch.float64, device=self.device)
        self.ntp_stoichiometry = torch.tensor([self.counts[base] for base in "AGCU"],
                                            dtype=torch.float64, device=self.device)

    def __getitem__(self, name):
        if name in self.parameters:
            return self.parameters[name]
        if name.startswith("nt_"):
            return self.counts[name[-1].upper()]
        if name == "AA" and self.construct:
            return sum(self.state[f"aa_{aa}"] for aa in AMINO_ACIDS)
        return self.state[name]

    def __setitem__(self, name, value):
        antimony_source(self.construct, {name: value} if name in PARAMETERS else None,
                        {name: value} if name not in PARAMETERS else None)
        if name in PARAMETERS:
            self.parameters[name] = value
        else:
            self.state[name] = value

    def resetAll(self):
        self.parameters = self.original_parameters.copy()
        self.state = self.original_initial.copy()

    def rhs(self, time, y):
        # Clamp flux inputs, not integrated states: validation sees negativity.
        return resource_rates(y, self.active_parameters, self.ntp_stoichiometry,
                              self.stoichiometry)

    def simulate_tensor(self, times, selections):
        t = torch.as_tensor(times, dtype=torch.float64, device=self.device)
        if t.ndim != 1 or len(t) < 2 or not torch.isfinite(t).all() or not (t[1:] > t[:-1]).all() or t[0] != 0:
            raise ValueError("Times must be finite, strictly increasing, and start at zero")
        y0 = torch.tensor([self.state[name] for name in self.species], dtype=torch.float64, device=self.device)
        self.active_parameters = torch.stack([
            torch.as_tensor(self.parameters[name], dtype=torch.float64, device=self.device)
            for name in PARAMETER_ORDER
        ])
        y = odeint(self.rhs, y0, t, method="dopri5",
                   rtol=self.integrator.relative_tolerance, atol=self.integrator.absolute_tolerance)
        if not torch.isfinite(y).all():
            raise RuntimeError("Nonfinite GPU trajectory")
        columns = []
        for name in selections:
            columns.append(t if name == "time" else y[:, 11:].sum(-1)
                           if name == "AA" and self.construct else y[:, self.species.index(name)])
        return torch.stack(columns, dim=1), y[-1]

    def simulate(self, start=0, end=240, points=241, selections=None, times=None):
        if times is None:
            if start != 0 or not math.isfinite(end) or end <= 0 or points < 2:
                raise ValueError("Start at zero with positive duration and at least two time points")
            times = np.linspace(start, end, points)
        selections = selections or ["time"] + self.species
        with torch.no_grad():
            values, final = self.simulate_tensor(times, selections)
        self.state = dict(zip(self.species, final.detach().cpu().tolist()))
        return Trajectory(values.detach().cpu().numpy(), selections)

    def getSBML(self):
        antimony.clearPreviousLoads()
        if antimony.loadAntimonyString(self.source) < 0:
            raise RuntimeError(antimony.getLastError())
        return antimony.getSBMLString("cfps")

    def execution_metadata(self):
        return dict(metadata(), solver="dopri5", flux_kernel="torch.compile", rtol=self.integrator.relative_tolerance,
                    atol=self.integrator.absolute_tolerance)


def load_model(construct=None, parameters=None, initial=None):
    return GPUModel(construct, parameters, initial)


def export_sbml(construct=None, parameters=None, initial=None):
    """Serialize the model without running a numerical backend."""
    model = antimony_source(construct, parameters, initial)
    antimony.clearPreviousLoads()
    if antimony.loadAntimonyString(model) < 0:
        raise RuntimeError(antimony.getLastError())
    return antimony.getSBMLString("cfps")


def simulate(minutes=240, points=241, construct=None, parameters=None, initial=None):
    if minutes <= 0 or points < 2:
        raise ValueError("Positive duration and at least two time points are required")
    rr = load_model(construct, parameters, initial)
    selections = ["time", "ATP", "GTP", "CTP", "UTP", "AA", "Fuel", "mRNA", "Protein", "GFP", "Machinery"]
    if construct is not None:
        selections += [f"aa_{aa}" for aa in AMINO_ACIDS]
    return rr, rr.simulate(0, minutes, points, selections=selections)
