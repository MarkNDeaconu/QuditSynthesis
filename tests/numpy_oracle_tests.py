"""Property tests for both backends against an external numpy oracle, with exact (coefficients, sde) canonical parity.

Run with:
    PYTHONPATH=src python tests/numpy_oracle_tests.py
"""

import math
import numpy as np
from quditsynthesis import datastructures as py
from quditsynthesis import datastructures_rust as rs

CASES = [
    (3, complex(0, math.sqrt(3))),
    (5, math.sqrt(5)),
    (7, complex(0, math.sqrt(7))),
    (8, math.sqrt(2)),
]


def direct_eval(coeffs, sde, p, loc):
    """Ground-truth complex value of Σ c_i ζ^i / loc^sde, independent of both backends."""
    zeta = np.exp(2j * np.pi / p)
    return sum(c * zeta ** i for i, c in enumerate(coeffs)) / complex(loc) ** sde


def rand_elem(ring, p, rng, max_c=5, min_sde=-2, max_sde=3):
    coeffs = [int(rng.integers(-max_c, max_c + 1)) for _ in range(p)]
    sde = int(rng.integers(min_sde, max_sde + 1))
    return py.cyclotomic_element(ring, coeffs, sde)


def to_rust(elem, ring_rust):
    return rs.cyclotomic_element_rust(ring_rust, elem.coefficients, elem.sde)


def canon(e):
    return (tuple(e.coefficients), e.sde)


def test_element_oracle(rng, ring, ring_rust, p, loc, n=200):
    """comp/conj/norm/add/mul against numpy ground truth, plus exact canonical parity."""
    for _ in range(n):
        raw_coeffs = [int(rng.integers(-5, 6)) for _ in range(p)]
        raw_sde = int(rng.integers(-2, 4))
        truth = direct_eval(raw_coeffs, raw_sde, p, loc)

        a = py.cyclotomic_element(ring, raw_coeffs, raw_sde)
        ar = rs.cyclotomic_element_rust(ring_rust, raw_coeffs, raw_sde)

        # Canonicalization preserves the value (both backends, vs oracle).
        assert np.isclose(a.comp(), truth), f"py comp {a.comp()} != truth {truth}"
        assert np.isclose(ar.comp(), truth), f"rs comp {ar.comp()} != truth {truth}"
        # Exact canonical parity between backends.
        assert canon(a) == canon(ar), f"canonical mismatch: {canon(a)} vs {canon(ar)}"

        # Conjugation against numpy (this catches sign bugs shared by both backends).
        assert np.isclose(a.conj().comp(), np.conj(truth)), f"py conj != np.conj at {canon(a)}"
        assert np.isclose(ar.conj().comp(), np.conj(truth)), f"rs conj != np.conj at {canon(ar)}"
        # Conjugation is an involution.
        assert canon(a.conj().conj()) == canon(a)

        # Norms are |value|² and never negative.
        assert np.isclose(a.norm(), abs(truth) ** 2), f"py norm {a.norm()} != |{truth}|²"
        assert np.isclose(ar.norm(), abs(truth) ** 2), f"rs norm {ar.norm()} != |{truth}|²"

        b = rand_elem(ring, p, rng)
        br = to_rust(b, ring_rust)
        tb = direct_eval(b.coefficients, b.sde, p, loc)

        s_py, s_rs = a + b, ar + br
        assert np.isclose(s_py.comp(), truth + tb), "py add != oracle"
        assert canon(s_py) == canon(s_rs), "add canonical mismatch"

        m_py, m_rs = a * b, ar * br
        assert np.isclose(m_py.comp(), truth * tb), "py mul != oracle"
        assert canon(m_py) == canon(m_rs), "mul canonical mismatch"

        k = int(rng.integers(-4, 5))
        sm_py, sm_rs = a * k, ar * k
        assert np.isclose(sm_py.comp(), truth * k), "py scalar mul != oracle"
        assert canon(sm_py) == canon(sm_rs), "scalar mul canonical mismatch"


def test_scalar_contract(ring, ring_rust, p, loc):
    """Non-integer scalars: exact ±|λ|^k absorption, TypeError otherwise — both backends."""
    coeffs = [1] + [0] * (p - 1)
    a = py.cyclotomic_element(ring, coeffs, 0)
    ar = rs.cyclotomic_element_rust(ring_rust, coeffs, 0)
    mag = abs(complex(loc))

    for k in (-2, 2):
        scalar = mag ** k
        truth = direct_eval(coeffs, 0, p, loc) * scalar
        r_py, r_rs = a * scalar, ar * scalar
        assert np.isclose(r_py.comp(), truth), f"py |λ|^{k} absorption wrong"
        assert canon(r_py) == canon(r_rs), f"|λ|^{k} canonical mismatch"
        assert all(isinstance(c, int) for c in r_py.coefficients), "float leaked into coefficients"

    if p in (5, 8):
        # Odd powers (e.g. the (1/√p)·operator idiom) are exact for real λ.
        r_py, r_rs = a * (1 / mag), ar * (1 / mag)
        assert (canon(r_py) == canon(r_rs) == (tuple(coeffs), 1)), "1/|λ| must raise sde by 1"
    else:
        # √p is not a ring element for p ≡ 3 (mod 4).
        for backend_elem in (a, ar):
            try:
                backend_elem * (1 / mag)
                assert False, "odd power of √p must raise TypeError"
            except TypeError:
                pass

    for backend_elem in (a, ar):
        try:
            backend_elem * 2.5
            assert False, "2.5 must raise TypeError"
        except TypeError:
            pass


def gate_set(mod, ring, p):
    """Known unitaries per dimension, built identically for either backend."""
    E = (lambda c, s=0: py.cyclotomic_element(ring, c, s)) if mod is py else \
        (lambda c, s=0: rs.cyclotomic_element_rust(ring, c, s))
    Op = mod.operator if mod is py else rs.operator_rust
    unit = lambda i: [1 if j == i else 0 for j in range(p)]

    if p == 8:
        # Qubit H = (1/√2)[[1,1],[1,−1]] and T = diag(1, ζ_8).
        one, minus = E(unit(0), 1), E([-1] + [0] * 7, 1)
        H = Op(2, 2, [[one, one], [one, minus]])
        T = Op(2, 2, [[E(unit(0)), E([0] * 8)], [E([0] * 8), E(unit(1))]])
        return [H, T]
    d = p
    # H = DFT/λ: entry (j,k) = ζ^{jk}/λ.
    H = Op(d, d, [[E(unit((j * k) % d), 1) for k in range(d)] for j in range(d)])
    # S = diag(1, ζ, 1, ..., 1); R = diag(1, ..., 1, −1); X = cyclic shift.
    n = E([0] * d)
    S = Op(d, d, [[E(unit(1)) if i == j == 1 else (E(unit(0)) if i == j else n) for j in range(d)] for i in range(d)])
    R = Op(d, d, [[(E([-1] + [0] * (d - 1)) if i == d - 1 else E(unit(0))) if i == j else n for j in range(d)] for i in range(d)])
    X = Op(d, d, [[E(unit(0)) if (i - j) % d == 1 else n for j in range(d)] for i in range(d)])
    return [H, S, R, X]


def test_known_gates_unitary(ring, ring_rust, p, loc):
    """unitary_check must accept genuine unitaries in BOTH backends (incl. odd sde)."""
    gates_py = gate_set(py, ring, p)
    gates_rs = gate_set(rs, ring_rust, p)
    for g_py, g_rs in zip(gates_py, gates_rs):
        u = np.asarray(g_py.comp())
        assert np.allclose(u @ u.conj().T, np.eye(u.shape[0]), atol=1e-8), "gate is not numerically unitary"
        assert g_py.unitary_check(), f"py unitary_check rejected a unitary (p={p})"
        assert g_rs.unitary_check(), f"rs unitary_check rejected a unitary (p={p})"
        assert np.allclose(u, np.asarray(g_rs.comp())), "backends evaluate the gate differently"


def test_tensor_is_kron(ring, ring_rust, p, rng):
    dim = 2 if p == 8 else p
    gates_py = gate_set(py, ring, p)
    gates_rs = gate_set(rs, ring_rust, p)
    A_py, B_py = gates_py[0], gates_py[-1]
    A_rs, B_rs = gates_rs[0], gates_rs[-1]
    T_py = A_py.tensor(B_py)
    T_rs = A_rs.tensor(B_rs)
    expected = np.kron(np.asarray(A_py.comp()), np.asarray(B_py.comp()))
    assert np.allclose(np.asarray(T_py.comp()), expected), "py tensor != np.kron"
    assert np.allclose(np.asarray(T_rs.comp()), expected), "rs tensor != np.kron"
    assert (T_rs.m, T_rs.n) == (dim * dim, dim * dim)
    # Tensor of unitaries is unitary — this failed under the outer-product bug.
    assert T_rs.unitary_check(), "tensor of unitaries must pass unitary_check"


def test_inner_product_and_state_norm(ring, ring_rust, p, loc):
    d = 2 if p == 8 else p
    gates_py = gate_set(py, ring, p)
    gates_rs = gate_set(rs, ring_rust, p)
    # First column of H is a normalized state (entries 1/λ each).
    col_py = [gates_py[0].matrix[i][0] for i in range(d)]
    col_rs = [rs.cyclotomic_element_rust(ring_rust, e.coefficients, e.sde) for e in col_py]
    st_py = py.state(d, col_py)
    st_rs = rs.state_rust(d, col_rs)
    vec = np.array([e.comp() for e in col_py])
    truth = np.vdot(vec, vec)
    assert np.isclose(st_py.norm(), truth) and np.isclose(st_py.norm(), 1.0), f"py state.norm {st_py.norm()} != 1"
    assert np.isclose(st_rs.norm(), truth) and np.isclose(st_rs.norm(), 1.0), f"rs state.norm {st_rs.norm()} != 1"
    # operator.__mul__ vector branch: ⟨u|v⟩ as a 1×1 operator, both backends.
    ip_py = st_py * st_py
    ip_rs = st_rs * st_rs
    assert np.isclose(ip_py.comp()[0][0], truth), "py inner product != vdot"
    assert np.isclose(ip_rs.comp()[0][0], truth), "rs inner product != vdot"
    _ = gates_rs


def test_disguised_zero(ring, ring_rust, p):
    """k·(1+ζ+…+ζ^{p-1}) must canonicalize to THE zero (sde 0) in both backends."""
    zero_py = py.cyclotomic_element(ring, [0] * p, 0)
    zero_rs = rs.cyclotomic_element_rust(ring_rust, [0] * p, 0)
    for k, sde in ((1, 5), (2, 1), (5, 3)):
        if p == 8:
            coeffs = [k, 0, 0, 0, k, 0, 0, 0]  # k·(1 + ζ⁴) = 0
        else:
            coeffs = [k] * p
        d_py = py.cyclotomic_element(ring, coeffs, sde)
        d_rs = rs.cyclotomic_element_rust(ring_rust, coeffs, sde)
        assert d_py == zero_py and hash(d_py) == hash(zero_py), f"py disguised zero {coeffs}/{sde}"
        assert d_rs == zero_rs and hash(d_rs) == hash(zero_rs), f"rs disguised zero {coeffs}/{sde}"
        assert len({d_py, zero_py}) == 1 and len({d_rs, zero_rs}) == 1


def test_group_ops(ring, ring_rust, p):
    """Orbits of real, non-commuting generators must agree EXACTLY across backends."""
    gates_py = gate_set(py, ring, p)
    gates_rs = gate_set(rs, ring_rust, p)
    d = 2 if p == 8 else p
    # H and X (or H and T for qubits): non-commuting, sde-changing.
    gens_py = [gates_py[0], gates_py[-1]]
    gens_rs = [gates_rs[0], gates_rs[-1]]
    gens_py[0].string, gens_py[1].string = "a", "b"
    gens_rs[0].string, gens_rs[1].string = "a", "b"

    depth = 4
    orb_py = ring.subgroup_bfs(gens_py, depth)
    orb_rs = ring_rust.subgroup_bfs_rust(gens_rs, depth)
    key_py = sorted(tuple(canon(e) for e in op.matrix.flatten()) for op in orb_py)
    key_rs = sorted(tuple(canon(e) for e in op.matrix.flatten()) for op in orb_rs)
    assert key_py == key_rs, f"orbit mismatch: {len(orb_py)} vs {len(orb_rs)} (p={p})"

    zero_py = py.cyclotomic_element(ring, [0] * p, 0)
    one_py = py.cyclotomic_element(ring, [1] + [0] * (p - 1), 0)
    zero_rs = rs.cyclotomic_element_rust(ring_rust, [0] * p, 0)
    one_rs = rs.cyclotomic_element_rust(ring_rust, [1] + [0] * (p - 1), 0)
    assert len(ring.torus(orb_py, zero_py)) == len(ring_rust.torus_rust(orb_rs, zero_rs)), "torus size"
    assert len(ring.permutation_subgroup(orb_py, one_py, zero_py)) == \
        len(ring_rust.permutation_subgroup_rust(orb_rs, one_rs, zero_rs)), "permutation size"

    # synth_search parity on a product with raised sde.
    target_py = gens_py[0] * gens_py[1] * gens_py[0]
    target_rs = gens_rs[0] * gens_rs[1] * gens_rs[0]
    res_py = target_py.synth_search(gens_py)
    res_rs = target_rs.synth_search_rust(gens_rs)
    if res_py is None:
        assert res_rs is None, "synth_search: py None, rs found one"
    else:
        assert res_rs is not None and res_py[1] == res_rs[1], "synth_search string mismatch"
    _ = d


def test_multiply_many(rng, ring, ring_rust, p, n=10):
    for _ in range(n):
        size = int(rng.integers(2, 4))
        length = int(rng.integers(2, 8))
        ops_py, ops_rs = [], []
        for _ in range(length):
            mat = [[rand_elem(ring, p, rng) for _ in range(size)] for _ in range(size)]
            ops_py.append(py.operator(size, size, mat))
            ops_rs.append(rs.operator_rust(size, size, [[to_rust(e, ring_rust) for e in row] for row in mat]))
        prod_py = ops_py[0]
        for op in ops_py[1:]:
            prod_py = prod_py * op
        prod_rs = rs.multiply_many_rust(ops_rs)
        assert np.allclose(np.asarray(prod_py.comp()), np.asarray(prod_rs.comp())), "multiply_many mismatch"
        assert sorted(prod_py.sde_profile().flatten()) == sorted(prod_rs.sde_profile().flatten())


def test_error_paths(ring_rust, p):
    """User mistakes raise catchable Python exceptions — never a process abort."""
    E = lambda c, s=0: rs.cyclotomic_element_rust(ring_rust, c, s)
    one = E([1] + [0] * (p - 1))

    try:
        _rust = rs._rust
        _rust.CyclotomicElementRust(p, [1] * (p + 1), 0)
        assert False, "wrong coefficient count must raise"
    except ValueError:
        pass

    if p != 5:
        other = rs.cyclotomic_ring_rust(5, math.sqrt(5))
        try:
            one + rs.cyclotomic_element_rust(other, [1, 0, 0, 0, 0], 0)
            assert False, "dim mismatch must raise"
        except ValueError:
            pass

    d = 2 if p == 8 else p
    v = rs.operator_rust(d, 1, [[one]] + [[E([0] * p)]] * (d - 1))
    sq = rs.operator_rust(d, d, [[one if i == j else E([0] * p) for j in range(d)] for i in range(d)])
    try:
        v * sq  # (d,1) x (d,d): invalid
        assert False, "shape mismatch must raise"
    except ValueError:
        pass

    # Overflow unwinds into a catchable exception; the interpreter survives.
    try:
        E([2 ** 61] + [0] * (p - 1)) * 8
        assert False, "overflow must raise"
    except BaseException as e:  # pyo3 PanicException derives from BaseException
        assert "overflow" in str(e).lower()

    try:
        sq.synthesize_rust([sq], target_sde=-10)
        assert False, "stalled synthesize must raise RuntimeError"
    except RuntimeError:
        pass


def main():
    rng = np.random.default_rng(0)
    for p, loc in CASES:
        ring = py.cyclotomic_ring(p, loc)
        ring_rust = rs.cyclotomic_ring_rust(p, loc)
        print(f"testing p={p}")
        test_element_oracle(rng, ring, ring_rust, p, loc)
        print("  element oracle OK (comp/conj/norm/add/mul vs numpy + exact parity)")
        test_scalar_contract(ring, ring_rust, p, loc)
        print("  scalar contract OK")
        test_known_gates_unitary(ring, ring_rust, p, loc)
        print("  known gates unitary OK")
        test_tensor_is_kron(ring, ring_rust, p, rng)
        print("  tensor == kron OK")
        test_inner_product_and_state_norm(ring, ring_rust, p, loc)
        print("  inner product / state norm OK")
        test_disguised_zero(ring, ring_rust, p)
        print("  disguised zero OK")
        test_group_ops(ring, ring_rust, p)
        print("  group ops OK (exact orbit parity)")
        test_multiply_many(rng, ring, ring_rust, p)
        print("  multiply_many OK")
        test_error_paths(ring_rust, p)
        print("  error paths OK (exceptions, not aborts)")
    print("ALL OK")


if __name__ == "__main__":
    main()
