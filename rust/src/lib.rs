mod element;
mod group;
mod operator;
mod ring;

use element::CyclotomicElement;
use operator::Operator;
use pyo3::exceptions::{PyIndexError, PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
use ring::RingDim;

fn ring_dim(d: usize) -> PyResult<RingDim> {
    RingDim::from_usize(d).ok_or_else(|| {
        PyValueError::new_err(format!("unsupported dimension {d} (expected 3, 5, 7 or 8)"))
    })
}

fn coeffs_from_list(dim: RingDim, list: &Bound<'_, PyList>) -> PyResult<[i64; 8]> {
    if list.len() != dim.value() {
        return Err(PyValueError::new_err(format!(
            "expected {} coefficients for dimension {}, got {}",
            dim.value(),
            dim.value(),
            list.len()
        )));
    }
    let mut coeffs = [0i64; 8];
    for (i, item) in list.iter().enumerate() {
        coeffs[i] = item
            .extract::<i64>()
            .map_err(|_| PyTypeError::new_err("coefficients must be integers"))?;
    }
    Ok(coeffs)
}

/// Canonical localization as (re, im): the Gauss sum g_p for odd primes, √2 for p=8.
fn localization_value(dim: RingDim) -> (f64, f64) {
    match dim {
        RingDim::D3 => (0.0, 3.0f64.sqrt()),
        RingDim::D5 => (5.0f64.sqrt(), 0.0),
        RingDim::D7 => (0.0, 7.0f64.sqrt()),
        RingDim::D8 => (2.0f64.sqrt(), 0.0),
    }
}

fn check_same_dim(a: RingDim, b: RingDim) -> PyResult<()> {
    if a != b {
        return Err(PyValueError::new_err(format!(
            "ring dimension mismatch: {} vs {}",
            a.value(),
            b.value()
        )));
    }
    Ok(())
}

#[pyclass(name = "CyclotomicElementRust")]
#[derive(Clone)]
struct PyCyclotomicElement {
    inner: CyclotomicElement,
}

#[pymethods]
impl PyCyclotomicElement {
    #[new]
    #[pyo3(signature = (dim, coeffs, sde=0))]
    fn new(dim: usize, coeffs: &Bound<'_, PyList>, sde: i32) -> PyResult<Self> {
        let dim = ring_dim(dim)?;
        let coeffs = coeffs_from_list(dim, coeffs)?;
        Ok(Self {
            inner: CyclotomicElement::new(dim, coeffs, sde),
        })
    }

    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim.value()
    }

    #[getter]
    fn coefficients(&self) -> Vec<i64> {
        self.inner.coeffs[0..self.inner.dimension()].to_vec()
    }

    #[getter]
    fn sde(&self) -> i32 {
        self.inner.sde
    }

    fn __add__(&self, other: &Bound<'_, PyCyclotomicElement>) -> PyResult<PyCyclotomicElement> {
        let other = other.borrow();
        check_same_dim(self.inner.dim, other.inner.dim)?;
        Ok(PyCyclotomicElement {
            inner: self.inner.add(&other.inner),
        })
    }

    fn __mul__<'py>(&self, py: Python<'py>, other: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(other) = other.downcast::<PyCyclotomicElement>() {
            let other = other.borrow();
            check_same_dim(self.inner.dim, other.inner.dim)?;
            let result = Bound::new(
                py,
                PyCyclotomicElement {
                    inner: self.inner.mul(&other.inner),
                },
            )?;
            return Ok(result.into_any());
        }
        if let Ok(scalar) = other.extract::<i64>() {
            let result = Bound::new(
                py,
                PyCyclotomicElement {
                    inner: self.inner.mul_scalar(scalar),
                },
            )?;
            return Ok(result.into_any());
        }
        Err(PyTypeError::new_err(
            "unsupported operand type for *: expected CyclotomicElementRust or int",
        ))
    }

    fn scale_localization(&self, sign: i64, k: i32) -> PyResult<PyCyclotomicElement> {
        self.inner
            .scale_localization(sign, k)
            .map(|inner| PyCyclotomicElement { inner })
            .map_err(PyValueError::new_err)
    }

    fn conj(&self) -> PyCyclotomicElement {
        PyCyclotomicElement {
            inner: self.inner.conj(),
        }
    }

    fn is_monomial(&self) -> bool {
        self.inner.is_monomial()
    }

    fn comp(&self) -> (f64, f64) {
        self.inner.comp(localization_value(self.inner.dim))
    }

    fn norm(&self) -> f64 {
        let c = self.inner.comp(localization_value(self.inner.dim));
        let cc = self.inner.conj().comp(localization_value(self.inner.dim));
        c.0 * cc.0 - c.1 * cc.1
    }

    fn __repr__(&self) -> String {
        format!(
            "CyclotomicElementRust(dim={}, coeffs={:?}, sde={})",
            self.inner.dim.value(),
            &self.inner.coeffs[0..self.inner.dimension()],
            self.inner.sde
        )
    }

    fn __eq__(&self, other: &Bound<'_, PyCyclotomicElement>) -> bool {
        self.inner == other.borrow().inner
    }

    fn __hash__(&self) -> u64 {
        self.inner.hash_value()
    }
}

#[pyclass(name = "OperatorRust")]
#[derive(Clone)]
struct PyOperator {
    inner: Operator,
}

fn operators_from_list(list: &Bound<'_, PyList>) -> PyResult<Vec<Operator>> {
    let mut ops = Vec::with_capacity(list.len());
    for item in list.iter() {
        let py_op = item
            .downcast::<PyOperator>()
            .map_err(|_| PyTypeError::new_err("expected a list of OperatorRust"))?;
        ops.push(py_op.borrow().inner.clone());
    }
    if let Some(first) = ops.first() {
        for op in &ops[1..] {
            check_same_dim(first.dim, op.dim)?;
        }
    }
    Ok(ops)
}

#[pymethods]
impl PyOperator {
    #[new]
    #[pyo3(signature = (dim, m, n, entries, gate_string=None))]
    fn new(
        dim: usize,
        m: usize,
        n: usize,
        entries: &Bound<'_, PyList>,
        gate_string: Option<String>,
    ) -> PyResult<Self> {
        let dim = ring_dim(dim)?;
        if entries.len() != m * n {
            return Err(PyValueError::new_err(format!(
                "expected {} entries for a {}x{} operator, got {}",
                m * n,
                m,
                n,
                entries.len()
            )));
        }
        let mut elems = Vec::with_capacity(m * n);
        for item in entries.iter() {
            let py_elem = item
                .downcast::<PyCyclotomicElement>()
                .map_err(|_| PyTypeError::new_err("entries must be CyclotomicElementRust"))?;
            let py_elem = py_elem.borrow();
            check_same_dim(dim, py_elem.inner.dim)?;
            elems.push(py_elem.inner.clone());
        }
        Ok(Self {
            inner: Operator::new(dim, m, n, elems, gate_string.unwrap_or_default()),
        })
    }

    #[getter]
    fn dim(&self) -> usize {
        self.inner.dim.value()
    }

    #[getter]
    fn shape(&self) -> (usize, usize) {
        (self.inner.m, self.inner.n)
    }

    #[getter]
    fn gate_string(&self) -> String {
        self.inner.gate_string.clone()
    }

    #[setter]
    fn set_gate_string(&mut self, s: String) {
        self.inner.gate_string = s;
    }

    #[getter]
    fn sde(&self) -> i32 {
        self.inner.sde()
    }

    /// The (0,1) entry's sde, or 0 for single-column operators — mirrors the
    /// reference's `elements[0][1]` with its IndexError-to-0 fallback.
    #[getter]
    fn sde2(&self) -> i32 {
        if self.inner.n > 1 {
            self.inner.get(0, 1).sde
        } else {
            0
        }
    }

    fn get(&self, row: usize, col: usize) -> PyResult<PyCyclotomicElement> {
        if row >= self.inner.m || col >= self.inner.n {
            return Err(PyIndexError::new_err("index out of bounds"));
        }
        Ok(PyCyclotomicElement {
            inner: self.inner.get(row, col).clone(),
        })
    }

    fn sde_sum(&self) -> i64 {
        self.inner.sde_sum()
    }

    fn sde_profile(&self) -> Vec<i32> {
        self.inner.sde_profile()
    }

    fn __mul__<'py>(&self, py: Python<'py>, other: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(other) = other.downcast::<PyOperator>() {
            let other = other.borrow();
            check_same_dim(self.inner.dim, other.inner.dim)?;
            // Two column vectors of equal length: Hermitian inner product,
            // mirroring the reference operator.__mul__.
            let result = if self.inner.n == 1 && other.inner.n == 1 && self.inner.m == other.inner.m
            {
                self.inner.inner_product(&other.inner)
            } else {
                if self.inner.n != other.inner.m {
                    return Err(PyValueError::new_err(format!(
                        "shape mismatch for operator product: ({}, {}) x ({}, {})",
                        self.inner.m, self.inner.n, other.inner.m, other.inner.n
                    )));
                }
                self.inner.matmul(&other.inner)
            };
            let py_result = Bound::new(py, PyOperator { inner: result })?;
            return Ok(py_result.into_any());
        }
        if let Ok(scalar) = other.extract::<i64>() {
            let py_result = Bound::new(
                py,
                PyOperator {
                    inner: self.inner.mul_scalar(scalar),
                },
            )?;
            return Ok(py_result.into_any());
        }
        Err(PyTypeError::new_err(
            "unsupported operand type for *: expected OperatorRust or int",
        ))
    }

    fn scale_localization(&self, sign: i64, k: i32) -> PyResult<PyOperator> {
        self.inner
            .scale_localization(sign, k)
            .map(|inner| PyOperator { inner })
            .map_err(PyValueError::new_err)
    }

    fn tensor(&self, other: &Bound<'_, PyOperator>) -> PyResult<PyOperator> {
        let other = other.borrow();
        check_same_dim(self.inner.dim, other.inner.dim)?;
        Ok(PyOperator {
            inner: self.inner.tensor(&other.inner),
        })
    }

    fn dag(&self) -> PyOperator {
        PyOperator {
            inner: self.inner.dag(),
        }
    }

    fn comp(&self) -> Vec<(f64, f64)> {
        self.inner.comp(localization_value(self.inner.dim))
    }

    #[pyo3(signature = (tol=None))]
    fn unitary_check(&self, tol: Option<f64>) -> bool {
        let tol = tol.unwrap_or(1e-8);
        self.inner
            .unitary_check(localization_value(self.inner.dim), tol)
    }

    fn is_diag(&self, null_elem: &Bound<'_, PyCyclotomicElement>) -> bool {
        self.inner.is_diag(&null_elem.borrow().inner)
    }

    fn is_permutation(
        &self,
        one_elem: &Bound<'_, PyCyclotomicElement>,
        null_elem: &Bound<'_, PyCyclotomicElement>,
    ) -> bool {
        self.inner
            .is_permutation(&one_elem.borrow().inner, &null_elem.borrow().inner)
    }

    fn monomial_check(&self) -> bool {
        self.inner.monomial_check()
    }

    fn __eq__(&self, other: &Bound<'_, PyOperator>) -> bool {
        self.inner == other.borrow().inner
    }

    fn __hash__(&self) -> u64 {
        self.inner.hash_value()
    }

    fn __repr__(&self) -> String {
        format!(
            "OperatorRust(dim={}, shape=({}, {}), sde_sum={})",
            self.inner.dim.value(),
            self.inner.m,
            self.inner.n,
            self.inner.sde_sum()
        )
    }

    fn synth_search(
        &self,
        py: Python<'_>,
        dropping_set: &Bound<'_, PyList>,
    ) -> PyResult<Option<(PyOperator, String)>> {
        let set = operators_from_list(dropping_set)?;
        if let Some(first) = set.first() {
            check_same_dim(self.inner.dim, first.dim)?;
        }
        let op = &self.inner;
        let result = py.detach(|| group::synth_search(op, &set));
        Ok(result.map(|(op, s)| (PyOperator { inner: op }, s)))
    }

    #[pyo3(signature = (dropping_set, target_sde=None))]
    fn synthesize(
        &self,
        py: Python<'_>,
        dropping_set: &Bound<'_, PyList>,
        target_sde: Option<i32>,
    ) -> PyResult<String> {
        let set = operators_from_list(dropping_set)?;
        if let Some(first) = set.first() {
            check_same_dim(self.inner.dim, first.dim)?;
        }
        let target = target_sde.unwrap_or(1);
        let op = &self.inner;
        py.detach(|| group::synthesize(op, &set, target))
            .ok_or_else(|| {
                PyRuntimeError::new_err(format!(
                    "synthesize: no dropping gate reduces SDE below {target}"
                ))
            })
    }
}

#[pyfunction]
#[pyo3(signature = (generators, depth=10))]
fn subgroup_bfs_rust(
    py: Python<'_>,
    generators: &Bound<'_, PyList>,
    depth: usize,
) -> PyResult<Vec<PyOperator>> {
    let gens = operators_from_list(generators)?;
    let orbit = py.detach(|| group::subgroup_bfs(&gens, depth));
    Ok(orbit.into_iter().map(|op| PyOperator { inner: op }).collect())
}

#[pyfunction]
fn torus_rust(
    py: Python<'_>,
    subgroup: &Bound<'_, PyList>,
    null_elem: &Bound<'_, PyCyclotomicElement>,
) -> PyResult<Vec<PyOperator>> {
    let ops = operators_from_list(subgroup)?;
    let null = null_elem.borrow().inner.clone();
    let result = py.detach(|| group::torus(&ops, &null));
    Ok(result.into_iter().map(|op| PyOperator { inner: op }).collect())
}

#[pyfunction]
fn permutation_subgroup_rust(
    py: Python<'_>,
    subgroup: &Bound<'_, PyList>,
    one_elem: &Bound<'_, PyCyclotomicElement>,
    null_elem: &Bound<'_, PyCyclotomicElement>,
) -> PyResult<Vec<PyOperator>> {
    let ops = operators_from_list(subgroup)?;
    let one = one_elem.borrow().inner.clone();
    let null = null_elem.borrow().inner.clone();
    let result = py.detach(|| group::permutation_subgroup(&ops, &one, &null));
    Ok(result.into_iter().map(|op| PyOperator { inner: op }).collect())
}

#[pyfunction]
#[pyo3(signature = (g, h, right=true))]
fn quotient_rust(
    py: Python<'_>,
    g: &Bound<'_, PyList>,
    h: &Bound<'_, PyList>,
    right: bool,
) -> PyResult<Vec<PyOperator>> {
    let g_ops = operators_from_list(g)?;
    let h_ops = operators_from_list(h)?;
    if let (Some(a), Some(b)) = (g_ops.first(), h_ops.first()) {
        check_same_dim(a.dim, b.dim)?;
    }
    let result = py.detach(|| group::quotient(&g_ops, &h_ops, right));
    Ok(result.into_iter().map(|op| PyOperator { inner: op }).collect())
}

#[pyfunction]
fn multiply_many_rust(py: Python<'_>, operators: &Bound<'_, PyList>) -> PyResult<PyOperator> {
    let ops = operators_from_list(operators)?;
    for pair in ops.windows(2) {
        if pair[0].n != pair[1].m {
            return Err(PyValueError::new_err(format!(
                "shape mismatch in operator chain: ({}, {}) x ({}, {})",
                pair[0].m, pair[0].n, pair[1].m, pair[1].n
            )));
        }
    }
    py.detach(|| group::multiply_many(&ops))
        .map(|op| PyOperator { inner: op })
        .ok_or_else(|| PyValueError::new_err("empty operator list"))
}

/// Left-multiplication walk over `generators` driven by `indices`; returns every
/// prefix product [g[i0], g[i1]*g[i0], ...]. Batches random-walk sampling into
/// one FFI call — the caller picks the random indices cheaply in Python.
#[pyfunction]
fn multiply_selected_rust(
    py: Python<'_>,
    generators: &Bound<'_, PyList>,
    indices: Vec<usize>,
) -> PyResult<Vec<PyOperator>> {
    let gens = operators_from_list(generators)?;
    let n = gens.len();
    let result = py.detach(|| group::multiply_selected_chain(&gens, &indices))
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "empty index list or index out of range for {n} generators"
            ))
        })?;
    Ok(result.into_iter().map(|op| PyOperator { inner: op }).collect())
}

/// k-fold Kronecker power of an operator, in one FFI call.
#[pyfunction]
fn tensor_power_rust(
    py: Python<'_>,
    op: &Bound<'_, PyOperator>,
    power: usize,
) -> PyResult<PyOperator> {
    let inner = op.borrow().inner.clone();
    py.detach(|| group::tensor_power(&inner, power))
        .map(|op| PyOperator { inner: op })
        .ok_or_else(|| PyValueError::new_err("power must be >= 1"))
}

/// Square-residue counts modulo the dimension (coefficients of the Gauss sum).
#[pyfunction]
fn gauss_sequence_rust(dim: usize) -> PyResult<Vec<i64>> {
    let dim = ring_dim(dim)?;
    Ok(ring::gauss_sequence(dim)[0..dim.value()].to_vec())
}

#[pymodule(name = "_rust")]
fn quditsynthesis_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCyclotomicElement>()?;
    m.add_class::<PyOperator>()?;
    m.add_function(wrap_pyfunction!(gauss_sequence_rust, m)?)?;
    m.add_function(wrap_pyfunction!(subgroup_bfs_rust, m)?)?;
    m.add_function(wrap_pyfunction!(torus_rust, m)?)?;
    m.add_function(wrap_pyfunction!(permutation_subgroup_rust, m)?)?;
    m.add_function(wrap_pyfunction!(quotient_rust, m)?)?;
    m.add_function(wrap_pyfunction!(multiply_many_rust, m)?)?;
    m.add_function(wrap_pyfunction!(multiply_selected_rust, m)?)?;
    m.add_function(wrap_pyfunction!(tensor_power_rust, m)?)?;
    Ok(())
}
