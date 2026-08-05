"""Rust-backed mirrors of the `datastructures.py` classes (same API, `_rust` suffix).

All arithmetic is delegated to the compiled `quditsynthesis._rust` extension (PyO3);
the wrappers hold no copied state — attribute reads are live views of the Rust objects.
"""

import math
import random
import numpy as np
from tabulate import tabulate

try:
    from quditsynthesis import _rust
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "The quditsynthesis._rust extension is not built. Build it with:\n"
        "  pip install maturin && maturin develop --release\n"
        "(run from the repository root; requires a Rust toolchain)"
    ) from e

from quditsynthesis.datastructures import superscript_map

# Reduction/comp conventions are fixed by the Gauss-sum localization (√2 for p=8); any
# other value would silently change the meaning of the sde, so ring construction validates it.
_CANONICAL_LOCALIZATION = {
    3: complex(0, math.sqrt(3)),
    5: complex(math.sqrt(5), 0),
    7: complex(0, math.sqrt(7)),
    8: complex(math.sqrt(2), 0),
}


def gauss_sequence_rust(p):
    """Rust-backed square-residue count (matches `gauss_sequence`)."""
    return _rust.gauss_sequence_rust(p)


def multiply_many_rust(operators):
    """Multiply a list of `operator_rust` matrices left-to-right; the whole chain runs in one FFI call."""
    if not operators:
        raise ValueError("multiply_many_rust requires at least one operator")
    inners = [op._inner for op in operators]
    return operator_rust._wrap(_rust.multiply_many_rust(inners), operators[0].ring)


def multiply_selected_rust(generators, indices):
    """Prefix products of the left-multiplication walk g[i_k]*...*g[i_1]*g[i_0]; the whole walk runs in one FFI call."""
    if not generators:
        raise ValueError("multiply_selected_rust requires at least one generator")
    inners = [g._inner for g in generators]
    ops = _rust.multiply_selected_rust(inners, list(indices))
    return [operator_rust._wrap(op, generators[0].ring) for op in ops]


class cyclotomic_ring_rust:
    """Ring descriptor mirroring `cyclotomic_ring`; the arithmetic lives in Rust."""

    def __init__(self, root_of_unity, localization) -> None:
        canonical = _CANONICAL_LOCALIZATION.get(root_of_unity)
        if canonical is None:
            raise ValueError(f"unsupported dimension {root_of_unity} (expected 3, 5, 7 or 8)")
        if not np.isclose(complex(localization), canonical):
            raise ValueError(
                f"localization {localization} does not match the ring convention "
                f"{canonical} for p={root_of_unity}; the Rust backend (and the sde "
                f"bookkeeping in both backends) is only meaningful for that value"
            )
        self.root_of_unity = root_of_unity
        self.localization = localization
        self.num_coefficient = root_of_unity

    def __eq__(self, value: object) -> bool:
        if isinstance(value, cyclotomic_ring_rust):
            return self.root_of_unity == value.root_of_unity and self.localization == value.localization
        return False

    def __hash__(self):
        return hash((self.root_of_unity, self.localization))

    def subgroup(self, generators, depth=10000):
        """Random-walk sampling of the generated subgroup; the walk's products run in one FFI call."""
        indices = [random.randrange(len(generators)) for _ in range(depth + 1)]
        # the reference excludes the initial generator pick (first prefix)
        return list(set(multiply_selected_rust(generators, indices)[1:]))

    def from_orbit(self, generator_set, depth=100):
        indices = [random.randrange(len(generator_set)) for _ in range(depth + 1)]
        return multiply_selected_rust(generator_set, indices)[-1]

    def subgroup_bfs_rust(self, generators, depth=10):
        """BFS closure under multiplication; the whole frontier loop runs in one FFI call."""
        gens = [g._inner for g in generators]
        ops = _rust.subgroup_bfs_rust(gens, depth)
        return [operator_rust._wrap(op, self) for op in ops]

    def torus_rust(self, subgroup, null_element):
        """Diagonal subgroup (Rust backend)."""
        ops = [g._inner for g in subgroup]
        result = _rust.torus_rust(ops, null_element._inner)
        return [operator_rust._wrap(op, self) for op in result]

    def permutation_subgroup_rust(self, subgroup, one_element, null_element=None):
        """Permutation subgroup (Rust backend)."""
        if null_element is None:
            null_element = cyclotomic_element_rust(self, [0] * self.num_coefficient, 0)
        ops = [g._inner for g in subgroup]
        result = _rust.permutation_subgroup_rust(ops, one_element._inner, null_element._inner)
        return [operator_rust._wrap(op, self) for op in result]

    def quotient_rust(self, G, H, right=True):
        """Coset representatives (Rust backend)."""
        g_ops = [g._inner for g in G]
        h_ops = [h._inner for h in H]
        result = _rust.quotient_rust(g_ops, h_ops, right)
        return [operator_rust._wrap(op, self) for op in result]

    # Reference-name aliases, so `ring.subgroup_bfs(...)` works on either backend.
    subgroup_bfs = subgroup_bfs_rust
    torus = torus_rust
    permutation_subgroup = permutation_subgroup_rust
    quotient = quotient_rust


class cyclotomic_element_rust:
    """Rust-backed cyclotomic element. API matches `cyclotomic_element`."""

    def __init__(self, ring: cyclotomic_ring_rust, coefficients, sde=0) -> None:
        self.ring = ring
        self._inner = _rust.CyclotomicElementRust(ring.root_of_unity, list(coefficients), sde)

    @classmethod
    def _wrap(cls, ring, inner):
        obj = object.__new__(cls)
        obj.ring = ring
        obj._inner = inner
        return obj

    @property
    def coefficients(self):
        return self._inner.coefficients

    @property
    def sde(self):
        return self._inner.sde

    def __add__(self, value: object) -> object:
        return cyclotomic_element_rust._wrap(self.ring, self._inner + value._inner)

    def __mul__(self, value: object) -> object:
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if value == 0:
                return cyclotomic_element_rust(self.ring, [0] * self.ring.num_coefficient, 0)
            if isinstance(value, float) and not value.is_integer():
                # Same contract as the reference: non-integer scalars must be ±|λ|^k, absorbed exactly into the sde.
                k = math.log(abs(value), abs(self.ring.localization))
                if not math.isclose(k, round(k)):
                    raise TypeError(f"cannot multiply ring element by {value}: not an integer or a power of the localization")
                sign = 1 if value > 0 else -1
                try:
                    inner = self._inner.scale_localization(sign, round(k))
                except ValueError as e:
                    raise TypeError(f"cannot multiply ring element by {value}: {e}") from None
                return cyclotomic_element_rust._wrap(self.ring, inner)
            return cyclotomic_element_rust._wrap(self.ring, self._inner * int(value))
        return cyclotomic_element_rust._wrap(self.ring, self._inner * value._inner)

    def __rmul__(self, value):
        return self * value

    def power(self, value):
        result = self
        for _ in range(value - 1):
            result = result * self
        return result

    def conj(self):
        return cyclotomic_element_rust._wrap(self.ring, self._inner.conj())

    def comp(self):
        re, im = self._inner.comp()
        return complex(re, im)

    def norm(self):
        return self._inner.norm()

    def is_monomial(self):
        return self._inner.is_monomial()

    def __eq__(self, other):
        if isinstance(other, cyclotomic_element_rust):
            return self._inner == other._inner
        return False

    def __hash__(self):
        return self._inner.__hash__()

    def __repr__(self):
        # Same polynomial rendering as the reference class.
        poly_string = ''
        coeffs = self.coefficients
        for index in range(self.ring.num_coefficient):
            c = coeffs[index]
            if c == 0:
                continue
            term = str(abs(c)) if index == 0 else str(abs(c)) + "ζ" + superscript_map.get(str(index))
            if c < 0:
                poly_string += ('-' if poly_string == '' else ' - ') + term
            else:
                poly_string += (term if poly_string == '' else ' + ' + term)
        return poly_string


class operator_rust:
    """Rust-backed operator matrix. API matches `operator`."""

    def __init__(self, m, n, elements, gate_string="") -> None:
        if isinstance(elements, np.ndarray):
            elements = elements.tolist()

        rust_entries = []
        ring = None
        for row in elements:
            for elem in row:
                if not isinstance(elem, cyclotomic_element_rust):
                    raise TypeError(f"operator_rust entries must be cyclotomic_element_rust, got {type(elem).__name__}")
                rust_entries.append(elem._inner)
                ring = elem.ring
        if ring is None:
            raise ValueError("operator_rust requires at least one element")

        self._bind(_rust.OperatorRust(ring.root_of_unity, m, n, rust_entries, gate_string), ring)

    def _bind(self, inner, ring):
        self._inner = inner
        self.ring = ring
        self.m, self.n = inner.shape
        self.shape = (self.m, self.n)
        self.dim = ring.root_of_unity
        return self

    @classmethod
    def _wrap(cls, inner, ring):
        return object.__new__(cls)._bind(inner, ring)

    @property
    def sde(self):
        return self._inner.sde

    @property
    def string(self):
        return self._inner.gate_string

    @string.setter
    def string(self, s):
        self._inner.gate_string = s

    @property
    def matrix(self):
        """The entries as a numpy object array of `cyclotomic_element_rust`."""
        return np.array(
            [[cyclotomic_element_rust._wrap(self.ring, self._inner.get(i, j)) for j in range(self.n)]
             for i in range(self.m)],
            dtype=object,
        )

    def power(self, exponent):
        if exponent < 1:
            raise ValueError("power must be >= 1")
        return multiply_many_rust([self] * exponent)

    def tensor(self, oper):
        return operator_rust._wrap(self._inner.tensor(oper._inner), self.ring)

    def tensor_power(self, power):
        """Repeated Kronecker product; the whole loop runs in one FFI call."""
        return operator_rust._wrap(_rust.tensor_power_rust(self._inner, power), self.ring)

    def __mul__(self, value):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if isinstance(value, float) and not value.is_integer():
                # Mirror the element-level contract: absorb ±|λ|^k exactly.
                k = math.log(abs(value), abs(self.ring.localization))
                if not math.isclose(k, round(k)):
                    raise TypeError(f"cannot multiply operator by {value}: not an integer or a power of the localization")
                sign = 1 if value > 0 else -1
                try:
                    inner = self._inner.scale_localization(sign, round(k))
                except ValueError as e:
                    raise TypeError(f"cannot multiply operator by {value}: {e}") from None
                return operator_rust._wrap(inner, self.ring)
            return operator_rust._wrap(self._inner * int(value), self.ring)
        if isinstance(value, operator_rust):
            # Rust dispatches matmul vs column-vector inner product and concatenates gate strings.
            return operator_rust._wrap(self._inner * value._inner, self.ring)
        raise TypeError(f"cannot multiply operator_rust with {type(value).__name__}")

    def __rmul__(self, value):
        if isinstance(value, (int, float)):
            return self * value
        return NotImplemented

    def sde_profile(self):
        return np.array(self._inner.sde_profile()).reshape(self.m, self.n)

    def sde_sum(self):
        return self._inner.sde_sum()

    def comp(self):
        values = self._inner.comp()
        return np.array([complex(re, im) for (re, im) in values]).reshape(self.m, self.n)

    def unitary_check(self, tol=1e-8):
        return self._inner.unitary_check(tol)

    def monomial_check(self):
        return self._inner.monomial_check()

    def synth_search_rust(self, dropping_set):
        """First gate that lowers total SDE (Rust backend); None if none does."""
        ds = [g._inner for g in dropping_set]
        out = self._inner.synth_search(ds)
        if out is None:
            return None
        op_inner, s = out
        return operator_rust._wrap(op_inner, self.ring), s

    def synthesize_rust(self, dropping_set, target_sde=1):
        """Iterative synthesis until min entry SDE ≤ target_sde; the whole loop runs in one FFI call."""
        ds = [g._inner for g in dropping_set]
        return self._inner.synthesize(ds, target_sde)

    def is_diag(self, null_element):
        return self._inner.is_diag(null_element._inner)

    def is_permutation(self, one_element, null_element):
        return self._inner.is_permutation(one_element._inner, null_element._inner)

    # Reference-name aliases.
    synth_search = synth_search_rust
    synthesize = synthesize_rust

    def __lt__(self, other):
        return self.sde_sum() < other.sde_sum()

    def __gt__(self, other):
        return self.sde_sum() > other.sde_sum()

    def __eq__(self, other):
        if isinstance(other, operator_rust):
            return self._inner == other._inner
        return False

    def __hash__(self):
        return self._inner.__hash__()

    def __repr__(self):
        # Same tabulate rendering as the reference operator.
        mat = self.matrix
        rows = self.m
        placement = rows // 2 - 1
        scalars = []
        for i in range(rows):
            if i == placement:
                scalars.append('√' + str(round((self.ring.localization ** 2).real)) + '^(-' + str(self.sde) + ')')
            else:
                scalars.append('')
        headers = [''] + [f'Column {i}' for i in range(1, self.n + 1)]
        matrix_with_scalars = np.column_stack((scalars, mat))
        return tabulate(matrix_with_scalars, headers, tablefmt='fancy_grid')


class state_rust(operator_rust):
    """Rust-backed state vector. API matches `state`."""

    def __init__(self, d, unit_vector) -> None:
        rows = [[elem] for elem in unit_vector]
        super().__init__(d, 1, rows)

    def norm(self):
        # ⟨ψ|ψ⟩ via the exact Hermitian inner product computed in Rust.
        return (self * self).comp()[0][0]
