
# Qudit Quantum Synthesis with Clifford+R Gateset

## Overview

This project provides a comprehensive implementation of **qudit quantum synthesis**, focusing on operations with **Clifford+R gates** for qudits of dimension 3, 5, and 7. The package allows users to perform complex quantum gate operations, search for gate sequences, reduce operators' **Smallest Denominator Exponent (SDE)**, and much more.

### Key Features:
- Operations with **cyclotomic rings** in qudit spaces $ Z[\zeta_p, 1/\sqrt{±p}] $
- **Cyclotomic element** arithmetic (addition, multiplication)
- **Matrix operators** for qudits (Hadamard, Rotation, and Phase gates)
- **Random operator generation** from a gateset
- **Synthesis search** for Clifford operators
- Precomputed **Clifford operator lookup tables** for qudit dimensions 3 and 5

---

## Dependencies

This library is built with **Python 3** and requires the following dependencies:

```bash
pip install numpy tabulate matplotlib networkx
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

- **3-dits**:
  ```python
  e0 = cyclotomic_element(z3, [1, 0, 0])  # Result: 1
  e1 = cyclotomic_element(z3, [0, 1, 0])  # Result: ζ
  e2 = cyclotomic_element(z3, [0, 0, 1])  # Result: ζ^2
  element1 = cyclotomic_element(z3, [2, 3, 4], 4) # Result: 2+ 3ζ + 4ζ^2
  element2 = cyclotomic_element(z3, [2, -7, -14], 4) # Result: (2 - 7ζ - 14ζ^2) / (sqrt(-3)^4)
  ```

- **5-dits**:
  ```python
  e0 = cyclotomic_element(z5, [1, 0, 0, 0, 0])  # Result: 1
  e1 = cyclotomic_element(z5, [0, 1, 0, 0, 0])  # Result: ζ
  e2 = cyclotomic_element(z5, [0, 0, 1, 0, 0])  # Result: ζ^2
  element1 = cyclotomic_element(z7, [6, 1, 0, 0, 2], sde=0) # Result: 6 + ζ + 2ζ^4
  element2 = cyclotomic_element(z7, [1, 1, -1, 0, 0], sde=2) # Result: (1 + ζ - ζ^2) / (sqrt(5)^2)
  ```

**Example (7-dits)**:
```python
element = cyclotomic_element(z7, [2, 1, 3, 0, 0, 1, 2], sde=1) # Result: (2 + ζ + 3ζ^2 + ζ^5 + 2ζ^6)) / (sqrt(-7)^2)
```

### Operations with Cyclotomic Elements

Cyclotomic elements support basic arithmetic operations such as **addition** and **multiplication**.

**Addition Example (5-dits)**:
```python
e0 = cyclotomic_element(z5, [1,2,0,-3,0])
e1 = cyclotomic_element(z5, [0,1,5,-2,0])
result = e0 + e1  # Result: 1 + 3ζ + 5ζ^2 - ζ^3
```

**Multiplication Example (5-dits)**:
```python
e0 = cyclotomic_element(z5, [1,-1,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,2,0])
result = e0 * e1  # Result: ζ -ζ^2+ ζ^3 -2ζ^4
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

# Matrix Examples for Qudit Operators

## Qutrits (Dimension 3)

### 1. Hadamard Gate \( H \) (Qutrits)
```python
H_qutrit = operator(3, 3, [
    [e0, e0, e0],
    [e0, e1, e2],
    [e0, e2, e1]
])
# Result: 3x3 Hadamard matrix
```

### 2. R Gate \( R \) (Qutrits)
```python
R_qutrit = operator(3, 3, [
    [e0, n, n],
    [n, e0, n],
    [n, n, -e0]
])
# Result: 3x3 R matrix
```

### 3. Phase Gate \( S \) (Qutrits)
```python
S_qutrit = operator(3, 3, [
    [e0, n, n],
    [n, e1, n],
    [n, n, e0]
])
# Result: 3x3 Phase matrix
```

---

## Qu5its (Dimension 5)

### 1. Hadamard Gate \( H \) (Qu5its)
```python
H_qu5it = operator(5, 5, [
    [e0, e0, e0, e0, e0],
    [e0, e1, e2, e3, e4],
    [e0, e2, e4, e1, e3],
    [e0, e3, e1, e4, e2],
    [e0, e4, e3, e2, e1]
])
# Result: 5x5 Hadamard matrix
```

### 2. R Gate \( R \) (Qu5its)
```python
R_qu5it = operator(5, 5, [
    [e0, n, n, n, n],
    [n, e0, n, n, n],
    [n, n, e0, n, n],
    [n, n, n, e0, n],
    [n, n, n, n, -e0]
])
# Result: 5x5 R matrix
```

### 3. Phase Gate \( S \) (Qu5its)
```python
S_qu5it = operator(5, 5, [
    [e0, n, n, n, n],
    [n, e1, n, n, n],
    [n, n, e3, n, n],
    [n, n, n, e1, n],
    [n, n, n, n, e0]
])
# Result: 5x5 Phase matrix
```

### Operator Multiplication

Operators can be multiplied if they have compatible dimensions, following standard matrix multiplication rules.

**Example 1**:
```python
 = H_qu5it * R_qu5it * H_qu5it 

```
```bash
Ex1 = 
╒═════════╤═════════════════════╤═════════════════════╤═════════════════════╤═════════════════════╤═════════════════════╕
│         │ Column 1            │ Column 2            │ Column 3            │ Column 4            │ Column 5            │
╞═════════╪═════════════════════╪═════════════════════╪═════════════════════╪═════════════════════╪═════════════════════╡
│         │ 3                   │ 2 + 2ζ¹ + 2ζ² + 2ζ³ │ -2ζ³                │ -2ζ²                │ -2ζ¹                │
├─────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┤
│ √5^(-2) │ 2 + 2ζ¹ + 2ζ² + 2ζ³ │ -2ζ³                │ -2ζ²                │ -2ζ¹                │ 3                   │
├─────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┤
│         │ -2ζ³                │ -2ζ²                │ -2ζ¹                │ 3                   │ 2 + 2ζ¹ + 2ζ² + 2ζ³ │
├─────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┤
│         │ -2ζ²                │ -2ζ¹                │ 3                   │ 2 + 2ζ¹ + 2ζ² + 2ζ³ │ -2ζ³                │
├─────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┼─────────────────────┤
│         │ -2ζ¹                │ 3                   │ 2 + 2ζ¹ + 2ζ² + 2ζ³ │ -2ζ³                │ -2ζ²                │
╘═════════╧═════════════════════╧═════════════════════╧═════════════════════╧═════════════════════╧═════════════════════╛
```
---

**Example 2**:
```python
 Ex2 = H_qutrit * H_qutrit * S_qutrit * R_qutrit * S_qutrit * H_qutrit

```
```bash
Ex2 = 
╒══════════╤════════════╤════════════╤════════════╕
│          │ Column 1   │ Column 2   │ Column 3   │
╞══════════╪════════════╪════════════╪════════════╡
│ √-3^(-1) │ -1         │ -1         │ -1         │
├──────────┼────────────┼────────────┼────────────┤
│          │ 1          │ -1 - 1ζ¹   │ 1ζ¹        │
├──────────┼────────────┼────────────┼────────────┤
│          │ 1 + 1ζ¹    │ -1         │ -1ζ¹       │
╘══════════╧════════════╧════════════╧════════════╛
```
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
This generates a random matrix in the orbit of H,S,R for 5-dits

```bash
╒══════════╤══════════════════════════════╤══════════════════════════════════╤════════════════════════════════╤════════════════════════════════╤════════════════════════════════╕
│          │ Column 1                     │ Column 2                         │ Column 3                       │ Column 4                       │ Column 5                       │
╞══════════╪══════════════════════════════╪══════════════════════════════════╪════════════════════════════════╪════════════════════════════════╪════════════════════════════════╡
│          │ 722 - 647ζ¹ + 664ζ² + 768ζ³  │ 1423 + 236ζ¹ + 82ζ² + 126ζ³      │ 1096 + 255ζ¹ + 742ζ² + 774ζ³   │ 203 - 561ζ¹ - 665ζ² - 55ζ³     │ -159 - 563ζ¹ - 93ζ² + 1847ζ³   │
├──────────┼──────────────────────────────┼──────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼────────────────────────────────┤
│ √5^(-10) │ -235 + 440ζ¹ - 92ζ² - 1066ζ³ │ -238 + 262ζ¹ - 5ζ² + 688ζ³       │ 2148 + 529ζ¹ + 1330ζ² + 1240ζ³ │ -460 + 750ζ¹ + 1299ζ² + 408ζ³  │ -1465 - 1406ζ¹ - 662ζ² - 470ζ³ │
├──────────┼──────────────────────────────┼──────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼────────────────────────────────┤
│          │ -707 - 361ζ¹ - 437ζ² + 7ζ³   │ -747 - 531ζ¹ - 1415ζ² + 1255ζ³   │ -646 - 1317ζ¹ - 1304ζ² - 776ζ³ │ 1243 + 670ζ¹ + 442ζ² + 602ζ³   │ 812 - 831ζ¹ + 204ζ² + 332ζ³    │
├──────────┼──────────────────────────────┼──────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼────────────────────────────────┤
│          │ 918 + 344ζ¹ - 1090ζ² + 675ζ³ │ -1492 - 1549ζ¹ - 1362ζ² - 1450ζ³ │ 426 + 461ζ¹ + 1996ζ² + 1344ζ³  │ 924 - 192ζ¹ + 600ζ² + 165ζ³    │ 144 - 214ζ¹ + 281ζ² - 184ζ³    │
├──────────┼──────────────────────────────┼──────────────────────────────────┼────────────────────────────────┼────────────────────────────────┼────────────────────────────────┤
│          │ 632 - 996ζ¹ - 1165ζ² - 54ζ³  │ 24 + 272ζ¹ + 510ζ² + 181ζ³       │ -74 - 803ζ¹ - 694ζ² - 182ζ³    │ -1520 - 117ζ¹ - 1106ζ² + 940ζ³ │ -22 - 656ζ¹ + 590ζ² - 615ζ³    │
╘══════════╧══════════════════════════════╧══════════════════════════════════╧════════════════════════════════╧════════════════════════════════╧════════════════════════════════╛
```
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
print(f"Reduced operator: {new_oper}")
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

# `neighbors_mat` Function Documentation

This function calculates a set of neighboring operators by applying a set of edges to the input operator (`mat`), with optional "free" Clifford operators applied on the right side of the edges. It utilizes `synth_search` to determine which edge lowers the **Smallest Denominator Exponent (SDE)** and uses that optimized edge in the resulting set of neighbors.

## Parameters:
- **`mat`**: `operator`  
  The base operator that serves as the starting point for calculating neighbors.

- **`edges`**: `list[operator]`  
  A list of operators representing the edges that will be applied to `mat`.

- **`edgesandcliffords`**: `list[operator]`  
  A list of operators that have Clifford operators applied to them on the right, used for synthesis optimization.

## Returns:
- **`list[operator]`**  
  A list of neighboring operators obtained by multiplying `edges` with `mat`. The edge that lowers the **SDE** is replaced with the optimized operator from `synth_search`, if found.

## Method:
1. **Synthesis Search**: The function uses `synth_search` to determine the edge in `edgesandcliffords` that lowers the SDE when applied to `mat`.
   
2. **Neighbor Calculation**: If an edge that lowers the SDE is found, it is used for the corresponding neighbor. The other edges are directly multiplied by `mat`.

3. **Error Handling**: If `synth_search` fails, all edges are simply multiplied by `mat` without optimization.


# Animation Script Documentation

This script visualizes the nodes and edges generated by the `neighbors_mat` function. It uses **Matplotlib** for plotting and **NetworkX** for managing the graph structure. The animation involves interactive elements, including mouse hover and click detection to update and explore neighboring nodes based on the current node.

## Functions:

- **`get_prepend_string(center_node, outer_node)`**  
  Returns a string based on the comparison of the `center_node` and `outer_node`, used to label the visualization.

- **`update_visualization(central_node)`**  
  Updates the graph visualization for the given `central_node`. It computes neighboring nodes using `neighbors_mat` and arranges them spatially, with the central node in the middle. The node labels and edges are drawn based on their properties.

- **`on_click(event)`**  
  Handles user clicks to select nodes on the plot. When a node is clicked, the visualization updates to display its neighbors.

- **`smooth_hover(event)`**  
  Handles mouse hover events, updating node labels to show additional information when a node is hovered.

## Key Details:
- **Interactive Visualization**: The graph is interactive, allowing users to click on nodes to explore neighbors or hover over them to display additional details.
- **Dynamic Layout**: The layout adjusts dynamically based on the node connections, positioning neighboring nodes either above or below the central node.
- **Hover and Click Events**: Mouse events are used to interact with the graph, providing a smooth user experience for navigating between nodes.


Feel free to reach out with any questions or feedback!
