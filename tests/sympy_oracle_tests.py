"""Exact oracle: both backends vs SymPy cyclotomic arithmetic.

Values are integer polynomials in x (= ζ_p) mod Φ_p, with λ substituted as its
exact polynomial in ζ (Gauss sum Σ (k|p)ζ^k for odd p, ζ+ζ⁷ = √2 for p=8), so
equality is exact polynomial congruence. Run: pytest tests/sympy_oracle_tests.py
"""

import math
import random

import pytest
import sympy as sp

from quditsynthesis import datastructures as py
from quditsynthesis import datastructures_rust as rs

x = sp.Symbol("x")

P_LOCS = {
    3: complex(0, math.sqrt(3)),
    5: complex(math.sqrt(5)),
    7: complex(0, math.sqrt(7)),
    8: complex(math.sqrt(2)),
}


class Oracle:
    def __init__(self, p):
        self.p = p
        self.Phi = sp.cyclotomic_poly(p, x)
        if p == 8:
            self.g = x + x**7
        else:
            self.g = sum(
                (1 if pow(k, (p - 1) // 2, p) == 1 else -1) * x**k
                for k in range(1, p)
            )

    def poly(self, coeffs):
        return sum(int(c) * x**i for i, c in enumerate(coeffs))

    def reduce(self, expr):
        return sp.rem(sp.expand(expr), self.Phi, x)

    def coeffs_of(self, expr):
        expr = sp.expand(expr)
        return [int(expr.coeff(x, i)) for i in range(self.p)]

    def eq(self, c1, s1, c2, s2):
        # c1/λ^s1 == c2/λ^s2  ⟺  c1·λ^s2 == c2·λ^s1 (mod Φ_p)
        p1, p2 = self.poly(c1), self.poly(c2)
        if s2 >= 0:
            p1 *= self.g**s2
        else:
            p2 *= self.g**(-s2)
        if s1 >= 0:
            p2 *= self.g**s1
        else:
            p1 *= self.g**(-s1)
        return self.reduce(p1 - p2) == 0


def make(backend, p):
    if backend == "python":
        ring = py.cyclotomic_ring(p, P_LOCS[p])
        return (lambda c, s=0: py.cyclotomic_element(ring, c, s)), py.operator
    ring = rs.cyclotomic_ring_rust(p, P_LOCS[p])
    return (lambda c, s=0: rs.cyclotomic_element_rust(ring, c, s)), rs.operator_rust


def rand_raw(rng, p):
    return [rng.randint(-4, 4) for _ in range(p)], rng.randint(-2, 3)


def canonical(elem):
    return [int(c) for c in elem.coefficients], int(elem.sde)


@pytest.mark.parametrize("p", [3, 5, 7, 8])
@pytest.mark.parametrize("backend", ["python", "rust"])
def test_reduction_preserves_value(backend, p):
    oracle = Oracle(p)
    E, _ = make(backend, p)
    rng = random.Random(f"{p}-{backend}")
    for _ in range(10):
        raw_c, raw_s = rand_raw(rng, p)
        c, s = canonical(E(raw_c, raw_s))
        assert oracle.eq(c, s, raw_c, raw_s)


@pytest.mark.parametrize("p", [3, 5, 7, 8])
@pytest.mark.parametrize("backend", ["python", "rust"])
def test_conj(backend, p):
    # Conjugation sends ζ→ζ⁻¹ and λ→-λ for p≡3(mod 4), hence (-1)^sde.
    oracle = Oracle(p)
    E, _ = make(backend, p)
    rng = random.Random(f"{p}-{backend}-conj")
    for _ in range(10):
        e = E(*rand_raw(rng, p))
        c0, s0 = canonical(e)
        conj = oracle.poly(c0).subs(x, x ** (p - 1))
        if p % 4 == 3 and s0 % 2 != 0:
            conj = -conj
        exp_c = oracle.coeffs_of(oracle.reduce(conj))
        c, s = canonical(e.conj())
        assert oracle.eq(c, s, exp_c, s0)


@pytest.mark.parametrize("p", [3, 5, 7, 8])
@pytest.mark.parametrize("backend", ["python", "rust"])
def test_matmul(backend, p):
    # Exercises element add/mul/reduction and the operator kernel together.
    oracle = Oracle(p)
    E, Op = make(backend, p)
    rng = random.Random(f"{p}-{backend}-matmul")
    for _ in range(5):
        A = [[E(*rand_raw(rng, p)) for _ in range(2)] for _ in range(2)]
        B = [[E(*rand_raw(rng, p)) for _ in range(2)] for _ in range(2)]
        C = Op(2, 2, A) * Op(2, 2, B)
        for i in range(2):
            for j in range(2):
                terms = [(A[i][k], B[k][j]) for k in range(2)]
                S = max(a.sde + b.sde for a, b in terms)
                total = sum(
                    oracle.poly(a.coefficients) * oracle.poly(b.coefficients)
                    * oracle.g ** (S - a.sde - b.sde)
                    for a, b in terms
                )
                c, s = canonical(C.matrix[i][j])
                assert oracle.eq(c, s, oracle.coeffs_of(oracle.reduce(total)), S)
