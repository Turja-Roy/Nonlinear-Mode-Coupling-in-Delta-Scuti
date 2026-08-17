import itertools

import numpy as np
import pytest

from coupling.angular import satisfies_selection_rules
from coupling.triplets import DETUNING_CUT_DIMLESS, enumerate_triplets

SMALL = dict(l_max=2, n_max=8)


@pytest.fixture(scope="module")
def cut(bg):
    return DETUNING_CUT_DIMLESS * bg.omega_dyn


def _canonical(t):
    return (t.sum_mode, tuple(sorted(t.pair)))


def _restrict(efs):
    return {
        k: v
        for k, v in efs.items()
        if k[0] <= SMALL["l_max"] and abs(k[1]) <= SMALL["n_max"]
    }


def _brute_force(efs, cut):
    out = set()
    for triple in itertools.combinations_with_replacement(sorted(efs), 3):
        for i in range(3):
            sum_mode = triple[i]
            pair = tuple(sorted(triple[:i] + triple[i + 1 :]))
            ls = (sum_mode[0], pair[0][0], pair[1][0])
            if not satisfies_selection_rules(*ls):
                continue
            delta = -efs[sum_mode].omega + efs[pair[0]].omega + efs[pair[1]].omega
            if abs(delta) < cut:
                out.add((sum_mode, pair))
    return out


def test_matches_brute_force(efs, cut):
    small = _restrict(efs)
    fast = {_canonical(t) for t in enumerate_triplets(small, cut, l_max=SMALL["l_max"])}
    assert fast == _brute_force(small, cut)


def test_no_duplicates(efs, cut):
    triplets = enumerate_triplets(efs, cut)
    assert len({_canonical(t) for t in triplets}) == len(triplets)


def test_all_pass_the_filters(efs, cut):
    for t in enumerate_triplets(efs, cut):
        assert satisfies_selection_rules(*t.ls)
        assert abs(t.delta) < cut
        assert t.delta == pytest.approx(sum(t.omega), rel=1e-12)
        assert t.omega[0] < 0 and t.omega[1] > 0 and t.omega[2] > 0


def test_all_eight_l_combinations_appear(efs, cut):
    found = {tuple(sorted(t.ls)) for t in enumerate_triplets(efs, cut)}
    assert len(found) == 8


def test_same_sign_triples_never_resonate(efs, cut):
    """Justifies enumerating only the one-minus sign class."""
    w = np.sort([ef.omega for ef in efs.values()])
    assert 3 * w[0] > cut


def test_m_combinations_are_nonzero(efs, cut):
    from coupling.angular import T

    for t in enumerate_triplets(efs, cut)[::200]:
        combos = t.m_combinations()
        assert combos
        for ms in combos:
            assert T(*t.ls, *ms) != 0.0
            assert sum(ms) == 0


def test_modes_carry_the_signs(efs, cut):
    t = enumerate_triplets(efs, cut)[0]
    modes = t.modes(efs, t.m_combinations()[0])
    assert modes[0].s == -1 and modes[0].omega < 0
    assert all(mo.s == +1 and mo.omega > 0 for mo in modes[1:])
    assert sum(mo.omega for mo in modes) == pytest.approx(t.delta, rel=1e-12)
