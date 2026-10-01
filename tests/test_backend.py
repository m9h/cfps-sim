import math
import pytest
import torch

from cfps_sim.backend import cuda_device
from cfps_sim.model import load_model, simulate


def test_missing_cuda_refuses_execution(monkeypatch):
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    with pytest.raises(RuntimeError, match="requires a CUDA GPU"):
        cuda_device()
    with pytest.raises(RuntimeError, match="requires a CUDA GPU"):
        simulate()


def test_invalid_inputs_are_rejected_before_device_allocation():
    with pytest.raises(ValueError, match="absent"):
        load_model(initial={"magnesium": 10})
    with pytest.raises(ValueError, match="positive"):
        load_model(parameters={"K_nt": 0})


@pytest.mark.gpu
def test_gpu_solver_refinement_and_location():
    model, trajectory = simulate()
    reference = model["GFP"]
    model.resetAll()
    model.integrator.relative_tolerance = 1e-12
    model.integrator.absolute_tolerance = 1e-16
    times = [0., 60., 120., 240.]
    values, final = model.simulate_tensor(times, ["time", "GFP"])
    assert values.is_cuda and final.is_cuda
    assert values.dtype == torch.float64
    assert float(values[-1, 1]) == pytest.approx(reference, rel=1e-5)


@pytest.mark.gpu
@pytest.mark.parametrize("function", ["demo", "cost_demo"])
def test_gpu_acquisition(function):
    pytest.importorskip("botorch")
    from cfps_sim import bo
    result = getattr(bo, function)()
    assert result["device"] == "cuda"
    assert result["dtype"] == "float64"
    assert all(0 <= value <= 1 for value in result["candidate_normalized"])
    if function == "demo":
        assert math.isfinite(result["predicted_mean"])
        assert result["predicted_std"] > 0
    else:
        assert math.isfinite(result["predicted_median_usd_per_g"])
        lower, upper = result["conditional_latent_90pct_usd_per_g"]
        assert 0 < lower <= upper
