"""Reduced, uncalibrated resource model. Concentrations: mM; time: minutes."""
import antimony
import roadrunner
import math
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


def load_model(construct=None, parameters=None, initial=None):
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
    antimony.clearPreviousLoads()
    if antimony.loadAntimonyString(model) < 0:
        raise RuntimeError(antimony.getLastError())
    rr = roadrunner.RoadRunner(antimony.getSBMLString("cfps"))
    rr.integrator.absolute_tolerance = 1e-12
    rr.integrator.relative_tolerance = 1e-8
    return rr


def simulate(minutes=240, points=241, construct=None, parameters=None, initial=None):
    if minutes <= 0 or points < 2:
        raise ValueError("Positive duration and at least two time points are required")
    rr = load_model(construct, parameters, initial)
    selections = ["time", "ATP", "GTP", "CTP", "UTP", "AA", "Fuel", "mRNA", "Protein", "GFP", "Machinery"]
    if construct is not None:
        selections += [f"aa_{aa}" for aa in AMINO_ACIDS]
    return rr, rr.simulate(0, minutes, points, selections=selections)
