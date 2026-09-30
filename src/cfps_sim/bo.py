"""CPU float64 BoTorch demonstration on a synthetic discrete candidate set."""
def demo(seed=42):
    import torch
    from botorch.models import SingleTaskGP
    from botorch.models.transforms.outcome import Standardize
    from botorch.fit import fit_gpytorch_mll
    from botorch.acquisition.analytic import LogExpectedImprovement
    from gpytorch.mlls import ExactMarginalLogLikelihood

    torch.manual_seed(seed)
    torch.set_num_threads(1)
    x = torch.rand(12, 2, dtype=torch.double, device="cpu")
    y = (1 - (x[:, :1] - 0.65)**2 - 2*(x[:, 1:] - 0.35)**2)
    y += 0.01 * torch.randn_like(y)
    gp = SingleTaskGP(x, y, outcome_transform=Standardize(m=1))
    fit_gpytorch_mll(ExactMarginalLogLikelihood(gp.likelihood, gp))
    grid = torch.linspace(0, 1, 21, dtype=torch.double, device="cpu")
    candidates = torch.cartesian_prod(grid, grid)
    with torch.no_grad():
        scores = LogExpectedImprovement(gp, best_f=y.max())(candidates.unsqueeze(1))
        best = candidates[scores.argmax()]
        posterior = gp.posterior(best.reshape(1, 2))
    return {"device": "cpu", "dtype": "float64", "seed": seed,
            "candidate_normalized": best.tolist(),
            "predicted_mean": posterior.mean.item(),
            "predicted_std": posterior.variance.sqrt().item(),
            "torch_version": torch.__version__, "mps_available": torch.backends.mps.is_available(),
            "note": "Synthetic objective; normalized grid is not a competition recipe."}
