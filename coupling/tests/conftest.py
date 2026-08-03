import pathlib

import pytest

from coupling.modes import load_model

MODELS = pathlib.Path(__file__).resolve().parents[2] / "models"
MODEL = MODELS / "dsct_M2.0"


def _load(path):
    if not any(path.glob("gyre/detail/detail.*.h5")):
        pytest.skip(f"{path.name} not built")
    return load_model(path)


@pytest.fixture(scope="session")
def model():
    return _load(MODEL)


@pytest.fixture(scope="session")
def bg(model):
    return model[0]


@pytest.fixture(scope="session")
def efs(model):
    return model[1]


@pytest.fixture(scope="session")
def poly3():
    """n = 3 polytrope: independent structure, no MESA input, Gamma_1 constant."""
    return _load(MODELS / "poly_n3")


@pytest.fixture(scope="session")
def poly0():
    """n = 0: constant density, and the f-mode is analytic."""
    out = _load(MODELS / "poly_n0")
    if (2, 0) not in out[1]:
        pytest.skip("poly_n0 incomplete: no l = 2 f-mode, rerun gyre_ad.in")
    return out
