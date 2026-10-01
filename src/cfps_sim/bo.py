"""CUDA float64 BoTorch demonstration on a synthetic discrete candidate set."""
from .backend import cuda_device, metadata
def demo(seed=42):
    import torch
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.optim.fit import fit_gpytorch_mll_torch
    from botorch.acquisition.analytic import LogExpectedImprovement
    from gpytorch.mlls import ExactMarginalLogLikelihood

    torch.manual_seed(seed)
    device = cuda_device()
    x = torch.rand(12, 2, dtype=torch.double, device=device)
    y = (1 - (x[:, :1] - 0.65)**2 - 2*(x[:, 1:] - 0.35)**2)
    y += 0.01 * torch.randn_like(y)
    gp = SingleTaskGP(x, y, outcome_transform=Standardize(m=1))
    fit_gpytorch_mll_torch(ExactMarginalLogLikelihood(gp.likelihood, gp), step_limit=2000)
    gp.eval()
    gp.likelihood.eval()
    grid = torch.linspace(0, 1, 21, dtype=torch.double, device=device)
    candidates = torch.cartesian_prod(grid, grid)
    with torch.no_grad():
        scores = LogExpectedImprovement(gp, best_f=y.max())(candidates.unsqueeze(1))
        best = candidates[scores.argmax()]
        posterior = gp.posterior(best.reshape(1, 2))
    return {**metadata(), "seed": seed,
            "candidate_normalized": best.tolist(),
            "predicted_mean": posterior.mean.item(),
            "predicted_std": posterior.variance.sqrt().item(),
            "note": "Synthetic objective; normalized grid is not a competition recipe."}


def cost_demo(seed=42):
    """Select on log(product mass / known reagent cost), not yield alone.

    The yield surface, volume-dependent stock-cost surrogate, and observations
    are illustrative. This does not connect recipes to calibrated predictions.
    """
    import math
    import torch
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.optim.fit import fit_gpytorch_mll_torch
    from botorch.acquisition.logei import qLogNoisyExpectedImprovement
    from botorch.acquisition.objective import GenericMCObjective
    from botorch.sampling.normal import SobolQMCNormalSampler
    from gpytorch.mlls import ExactMarginalLogLikelihood

    torch.manual_seed(seed)
    device = cuda_device()
    volume_l = 20e-6
    grid = torch.linspace(0, 1, 21, dtype=torch.double, device=device)
    candidates = torch.cartesian_prod(grid, grid)
    indices = torch.randperm(len(candidates), device=device)[:12]
    x = candidates[indices]

    def costs(points):
        return 0.024 + 0.002*points[..., 0] + 0.1*points[..., 1]**2

    def true_yield(points):
        return 0.1 + 0.9*torch.exp(-3*(points[..., 0]-0.65)**2 - 8*(points[..., 1]-0.7)**2)

    log_y = true_yield(x).log().unsqueeze(-1) + 0.08*torch.randn(12, 1, dtype=torch.double, device=device)
    gp = SingleTaskGP(x, log_y, outcome_transform=Standardize(m=1))
    fit_gpytorch_mll_torch(ExactMarginalLogLikelihood(gp.likelihood, gp), step_limit=2000)
    gp.eval()
    gp.likelihood.eval()
    # Maximizing this utility is equivalent pointwise to minimizing USD/g.
    # Subtract known cost for both baseline and candidate posterior samples.
    objective = GenericMCObjective(lambda samples, X: samples[..., 0] + math.log(volume_l) - costs(X).log())
    acquisition = qLogNoisyExpectedImprovement(
        gp, X_baseline=x, objective=objective,
        sampler=SobolQMCNormalSampler(torch.Size([128]), seed=seed), prune_baseline=False,
    )
    remaining = torch.ones(len(candidates), dtype=torch.bool, device=device)
    remaining[indices] = False
    choices = candidates[remaining]
    with torch.no_grad():
        scores = torch.cat([acquisition(chunk.unsqueeze(1)) for chunk in choices.split(64)])
        selected = choices[scores.argmax()]
        posterior = gp.posterior(selected.reshape(1, 2))
        median_titer = posterior.mean.exp().item()
        mean_log = posterior.mean.item()
        sd_log = posterior.variance.sqrt().item()
        cost = costs(selected).item()
        true_titers, true_costs = true_yield(candidates), costs(candidates)
        best_yield = int(true_titers.argmax())
        best_economics = int((true_titers*volume_l/true_costs).argmax())

    def oracle(index):
        return {"candidate_normalized": candidates[index].tolist(), "titer_g_l": true_titers[index].item(),
                "cost_usd": true_costs[index].item(),
                "usd_per_g": (true_costs[index]/(true_titers[index]*volume_l)).item()}

    return {"data_kind": "synthetic", "seed": seed, **metadata(),
            "objective": "maximize log(g/USD), equivalent pointwise to minimize USD/g",
            "candidate_normalized": selected.tolist(), "reaction_cost_usd": cost,
            "predicted_median_titer_g_l": median_titer,
            "predicted_median_usd_per_g": cost/(median_titer*volume_l),
            "conditional_latent_90pct_usd_per_g": [cost/(volume_l*math.exp(mean_log+1.6448536*sd_log)),
                                                      cost/(volume_l*math.exp(mean_log-1.6448536*sd_log))],
            "synthetic_oracle_max_yield": oracle(best_yield),
            "synthetic_oracle_min_usd_per_g": oracle(best_economics),
            "note": "Toy costs include a fixed reagent floor. Oracle values explain the example and never enter GP training or acquisition. Intervals are conditional GP latent intervals, not calibrated biological or observation-error intervals. No recipe mapping, minimum-yield constraint, or experimental budget policy is implemented."}
