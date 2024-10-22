# Qudit Quantum Synthesis with Clifford+R Gateset

## Overview

This project provides a comprehensive implementation of **qudit quantum synthesis**, focusing on operations with **Clifford+R gates** for qudits of dimension 3, 5, and 7. The package allows users to perform complex quantum gate operations, search for gate sequences, reduce operators' **Smallest Denominator Exponent (SDE)**, and much more.

### Key Features:
- Operations with **cyclotomic rings** in qudit spaces (\( Z[\zeta_p, 1/\sqrt{±p}] \))
- **Cyclotomic element** arithmetic (addition, multiplication)
- **Matrix operators** for qudits (Hadamard, Rotation, and Phase gates)
- **Random operator generation** from a gateset
- **Synthesis search** for Clifford operators
- Precomputed **Clifford operator lookup tables** for qudit dimensions 3 and 5

---

## Dependencies

This library is built with **Python 3** and requires the following dependencies:

```bash
pip install numpy tabulate
```

- **NumPy**: For handling matrix operations.
- **Tabulate**: For displaying matrices in a readable format.
- **Pickle**: For loading precomputed Clifford operators.

Make sure to also import `pickle` in your code:

```python
import pickle
```

---

## Cyclotomic Rings and Elements

### Constructing Cyclotomic Rings

Cyclotomic rings \( Z[\zeta_p, 1/\sqrt{±p}] \) can be created using the `cyclotomic_ring` class. For example:

```python
z3 = cyclotomic_ring(3, complex(0, math.sqrt(3)))
z5 = cyclotomic_ring(5, math.sqrt(5))
z7 = cyclotomic_ring(7, complex(0, math.sqrt(7)))
```

### Constructing Cyclotomic Elements

Cyclotomic elements are created using the `cyclotomic_element` class. These elements represent polynomials in powers of \( \zeta \), and they can also have a scaling factor determined by the **Smallest Denominator Exponent (SDE)**.

**Syntax**:
```python
element = cyclotomic_element(ring, coefficients, sde=0)
```

- **ring**: The cyclotomic ring (e.g., `z3`, `z5`, `z7`).
- **coefficients**: List of coefficients for powers of \( \zeta \).
- **sde**: (Optional) Smallest Denominator Exponent, which determines the power on the denominator.

**Example (7-dits)**:
```python
element = cyclotomic_element(z7, [2, 1, 3, 0, 0, 1, 2], sde=1)
```
This represents the polynomial \( 2 + \zeta + 3\zeta^2 + \zeta^5 + 2\zeta^6 \), scaled by \( rac{1}{\sqrt{-7}} \).

### Operations with Cyclotomic Elements

Cyclotomic elements support basic arithmetic operations such as **addition** and **multiplication**.

**Addition Example (5-dits)**:
```python
e0 = cyclotomic_element(z5, [1,0,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,0,0])
result = e0 + e1  # Result: 1 + ζ
```

**Multiplication Example (5-dits)**:
```python
e0 = cyclotomic_element(z5, [1,0,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,0,0])
result = e0 * e1  # Result: ζ
```

---

## Operators for Qudits

### Matrix Operators

Operators in the form of matrices are created using the `operator` class. Each entry in the matrix is a `cyclotomic_element`, and the first two arguments define the number of rows and columns of the matrix.

**Syntax**:
```python
op = operator(rows, cols, matrix)
```

- **rows**: Number of rows.
- **cols**: Number of columns.
- **matrix**: List of lists containing `cyclotomic_element` objects.

**Example (5-dits Hadamard Gate)**:
```python
H = operator(5, 5, [[e0, e0, e0, e0, e0], [e0, e1, e2, e3, e4], ...])
```

### Operator Multiplication

Operators can be multiplied if they have compatible dimensions, following standard matrix multiplication rules.

**Example (5-dits \( H 	imes H \))**:
```python
H2 = H * H  # Performs matrix multiplication of H with itself
```
For the 5-dits Hadamard operator, multiplying \( H 	imes H \) results in the identity matrix \( I_5 \).

---

## Random Operator Generation: `from_orbit`

This function generates a random operator by multiplying elements from a gateset up to a given depth.

```python
def from_orbit(generator_set, depth=100):
    curr = generator_set[0]
    for i in range(depth):
        curr = random.choice(generator_set) * curr
    return curr
```

**Example (5-dits)**:
```python
from_orbit([H, R, S], depth=100)
```
This generates a random matrix with a scaling factor of \( rac{1}{\sqrt{5}^{10}} \).

---

## Synthesis Search: `synth_search`

The `synth_search` function attempts to reduce the **Smallest Denominator Exponent (SDE)** of an operator by left-multiplying it with elements from a specified set.

```python
def synth_search(oper, dropping_set):
    for option in dropping_set:
        new_oper = option * oper
        if new_oper.sde < oper.sde:
            return new_oper, option.string
```

**Example**:
```python
new_oper, reducer = synth_search(oper, [H, R, S])
print(f"SDE reduced using: {reducer}")
```

---

## Precomputed Clifford Operators

Lookup tables for Clifford operators in 3 and 5 dimensions are stored in **pickle files** for fast access:

```python
with open('cliffords5.pkl', 'rb') as f:
    cliffords = pickle.load(f)
```

- **`cliffords5.pkl`**: Contains all Clifford operators in 5 dimensions.
- **`cliffords3.pkl`**: Contains all Clifford operators in 3 dimensions.

---

### Important Note

**`cliffords7.pkl`**, which contains all Clifford operators in 7 dimensions, is **too large to be posted on GitHub**. Please **email me directly** if you need access to this file.

---

Feel free to reach out with any questions or feedback!

