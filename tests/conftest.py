import pytest
import torch


def pytest_collection_modifyitems(config, items):
    if not torch.cuda.is_available():
        for item in items:
            if "gpu" in item.keywords:
                item.add_marker(pytest.mark.skip(reason="Requires CUDA GPU"))
