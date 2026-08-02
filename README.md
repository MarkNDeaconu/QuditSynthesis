# QuditSynthesis
A software package for quantum computing related calculations and simulations over cyclotomic unitary groups.

Exact arithmetic in the rings \( Z[\zeta_p, 1/\sqrt{\pm p}] \) for qudit dimensions \( p = 3, 5, 7, 8 \) — no floating point until you ask for it. Two backends behind one API: a Rust extension (default, 10–400× faster) and a pure-Python reference.

## Getting Started

Requires Python 3 and a Rust toolchain.

```
git clone https://github.com/MarkNDeaconu/QuditSynthesis
cd QuditSynthesis
pip install .
```

## Overview

- **Cyclotomic ring arithmetic** — elements \( \sum c_i \zeta^i / \lambda^{sde} \) with exact canonical reduction by the Smallest Denominator Exponent (SDE).
- **Matrix operators** — qudit gates (Hadamard, Phase, R) with exact matmul, tensor products, conjugates, unitarity checks, and pretty printing.
- **Group algorithms** — `subgroup_bfs`, `torus` (diagonals), `permutation_subgroup`, `quotient` (coset representatives).
- **Synthesis** — `synth_search` and `synthesize`: reduce an operator's SDE by left-multiplication from a dropping set and recover the gate string.

## Example

```python
import math
from quditsynthesis import cyclotomic_ring, cyclotomic_element, operator

z5 = cyclotomic_ring(5, math.sqrt(5))
E = lambda c, s=0: cyclotomic_element(z5, c, s)
unit = lambda i: [1 if j == i else 0 for j in range(5)]

# Hadamard H = DFT/√5
H = operator(5, 5, [[E(unit((j * k) % 5), 1) for k in range(5)] for j in range(5)])

H.unitary_check()          # True, exact
product = H * H * H        # exact ring arithmetic in every entry

word = z5.from_orbit([H], depth=100)   # random operator from the orbit of H
word.synthesize([H], target_sde=1)     # gate string reducing SDE to ≤ 1
```

Printing an operator renders the exact matrix with its SDE prefactor.

## Backends

```python
import quditsynthesis as qs
qs.set_backend("python")          # or set QUDITSYNTHESIS_BACKEND=python before importing
```

Set the backend before importing names from the package. For long gate words use `multiply_many([A, B, C])` instead of `A*B*C` — the whole chain runs in a single Rust call. The backend modules can also be imported directly (`quditsynthesis.datastructures`, `quditsynthesis.datastructures_rust`).

## Testing

Correctness is checked against two independent oracles, each on both backends:

```
pytest tests/sympy_oracle_tests.py   # exact symbolic arithmetic vs SymPy
python tests/numpy_oracle_tests.py   # numerical oracle + backend parity
```

These run in CI on every push to main.
