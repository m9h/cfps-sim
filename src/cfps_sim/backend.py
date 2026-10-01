"""GPU execution policy shared by integration, calibration, and acquisition."""
import torch


def cuda_device():
    if not torch.cuda.is_available():
        raise RuntimeError(
            "CFPS numerical execution requires a CUDA GPU and a CUDA-enabled "
            "PyTorch installation. Run on a GPU host or the Modal GPU runner."
        )
    return torch.device("cuda")


def metadata():
    device = cuda_device()
    return {"backend": "torch", "device": str(device), "dtype": "float64",
            "gpu": torch.cuda.get_device_name(device), "torch_version": torch.__version__,
            "cuda_version": torch.version.cuda}
