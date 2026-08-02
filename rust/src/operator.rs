use crate::element::{conv_i128, conv_pair_i128, CyclotomicElement};
use crate::ring::{sde_alignment_multiplier, RingDim};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug)]
pub struct Operator {
    pub dim: RingDim,
    pub m: usize,
    pub n: usize,
    pub entries: Vec<CyclotomicElement>,
    pub gate_string: String,
}

impl Operator {
    pub fn new(dim: RingDim, m: usize, n: usize, entries: Vec<CyclotomicElement>, gate_string: String) -> Self {
        assert_eq!(entries.len(), m * n);
        Self {
            dim,
            m,
            n,
            entries,
            gate_string,
        }
    }

    #[inline]
    pub fn idx(&self, row: usize, col: usize) -> usize {
        row * self.n + col
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> &CyclotomicElement {
        &self.entries[self.idx(row, col)]
    }

    pub fn sde(&self) -> i32 {
        // First entry's SDE is the representative, matching the Python code.
        if self.entries.is_empty() {
            0
        } else {
            self.entries[0].sde
        }
    }

    pub fn sde_sum(&self) -> i64 {
        self.entries.iter().map(|e| e.sde as i64).sum()
    }

    pub fn sde_profile(&self) -> Vec<i32> {
        self.entries.iter().map(|e| e.sde).collect()
    }

    /// Each output entry accumulates its inner-product terms as raw convolutions
    /// aligned to a common sde and is reduced once — canonically identical to
    /// reduce-per-term (reduction is confluent), but without k−1 intermediate
    /// canonicalizations per entry.
    pub fn matmul(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        assert_eq!(self.n, other.m);
        let p = self.dim.value();
        let denom = sde_alignment_multiplier(self.dim);
        let (m, n, k) = (self.m, other.n, self.n);
        let mut entries = Vec::with_capacity(m * n);

        for i in 0..m {
            for j in 0..n {
                let mut target = i32::MIN;
                for l in 0..k {
                    let a = self.get(i, l);
                    let b = other.get(l, j);
                    if a.is_zero() || b.is_zero() {
                        continue;
                    }
                    target = target.max(a.sde + b.sde);
                }
                if target == i32::MIN {
                    entries.push(CyclotomicElement::new(self.dim, [0i64; 8], 0));
                    continue;
                }
                let mut acc = [0i128; 8];
                for l in 0..k {
                    let a = self.get(i, l);
                    let b = other.get(l, j);
                    if a.is_zero() || b.is_zero() {
                        continue;
                    }
                    let mut term = conv_pair_i128(a, b);
                    for _ in 0..(target - (a.sde + b.sde)) {
                        term = conv_i128(&term, denom, p);
                    }
                    for t in 0..p {
                        acc[t] += term[t];
                    }
                }
                entries.push(CyclotomicElement::from_i128(self.dim, acc, target));
            }
        }

        let mut res = Self::new(self.dim, m, n, entries, String::new());
        res.gate_string = format!("{}{}", self.gate_string, other.gate_string);
        res
    }

    /// Hermitian inner product of two column vectors: a 1×1 operator holding
    /// Σ_i conj(u_i)·v_i, mirroring the reference operator.__mul__ vector branch.
    pub fn inner_product(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        assert_eq!((self.n, other.n), (1, 1));
        assert_eq!(self.m, other.m);
        let mut total = CyclotomicElement::new(self.dim, [0i64; 8], 0);
        for i in 0..self.m {
            total = total.add(&self.get(i, 0).conj().mul(other.get(i, 0)));
        }
        Self::new(self.dim, 1, 1, vec![total], String::new())
    }

    /// Kronecker product: entry ((i1,i2),(j1,j2)) = A[i1][j1]·B[i2][j2].
    pub fn tensor(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        let m = self.m * other.m;
        let n = self.n * other.n;
        let zero = CyclotomicElement::new(self.dim, [0i64; 8], 0);
        let mut entries = vec![zero; m * n];
        for i1 in 0..self.m {
            for j1 in 0..self.n {
                let a = self.get(i1, j1);
                for i2 in 0..other.m {
                    for j2 in 0..other.n {
                        let b = other.get(i2, j2);
                        entries[(i1 * other.m + i2) * n + (j1 * other.n + j2)] = a.mul(b);
                    }
                }
            }
        }
        Self::new(self.dim, m, n, entries, String::new())
    }

    pub fn mul_scalar(&self, scalar: i64) -> Self {
        let entries = self
            .entries
            .iter()
            .map(|e| e.mul_scalar(scalar))
            .collect();
        Self::new(self.dim, self.m, self.n, entries, String::new())
    }

    pub fn scale_localization(&self, sign: i64, k: i32) -> Result<Self, &'static str> {
        let mut entries = Vec::with_capacity(self.m * self.n);
        for e in &self.entries {
            entries.push(e.scale_localization(sign, k)?);
        }
        Ok(Self::new(self.dim, self.m, self.n, entries, String::new()))
    }

    pub fn dag(&self) -> Self {
        let mut entries = Vec::with_capacity(self.m * self.n);
        for j in 0..self.n {
            for i in 0..self.m {
                entries.push(self.get(i, j).conj());
            }
        }
        Self::new(self.dim, self.n, self.m, entries, String::new())
    }

    pub fn comp(&self, localization: (f64, f64)) -> Vec<(f64, f64)> {
        self.entries.iter().map(|e| e.comp(localization)).collect()
    }

    /// ‖(U·U† − I)_{ij}‖ ≤ atol + rtol·|I_{ij}| with rtol = 1e-5, mirroring
    /// numpy.allclose. Conjugation is numerical — (U†)_{lj} = conj(U_{jl}) on
    /// the evaluated matrix — so no algebraic dag() and no allocation per step.
    pub fn unitary_check(&self, localization: (f64, f64), atol: f64) -> bool {
        if self.m != self.n {
            return false;
        }
        const RTOL: f64 = 1e-5;
        let n = self.n;
        let u = self.comp(localization);
        for i in 0..n {
            for j in 0..n {
                let mut re = 0.0;
                let mut im = 0.0;
                for l in 0..n {
                    let (ar, ai) = u[i * n + l];
                    let (br, bi) = u[j * n + l];
                    re += ar * br + ai * bi;
                    im += ai * br - ar * bi;
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                let dr = re - expected;
                if (dr * dr + im * im).sqrt() > atol + RTOL * expected {
                    return false;
                }
            }
        }
        true
    }

    pub fn is_diag(&self, null_elem: &CyclotomicElement) -> bool {
        for i in 0..self.m {
            for j in 0..self.n {
                if i != j && self.get(i, j) != null_elem {
                    return false;
                }
            }
        }
        true
    }

    /// Exactly one `one_elem` per row, all other entries `null_elem`. Columns are
    /// deliberately not checked — this mirrors the Python reference exactly, and
    /// is sufficient for unitary inputs.
    pub fn is_permutation(&self, one_elem: &CyclotomicElement, null_elem: &CyclotomicElement) -> bool {
        for i in 0..self.m {
            let mut one_count = 0;
            for j in 0..self.n {
                let e = self.get(i, j);
                if e == one_elem {
                    one_count += 1;
                } else if e != null_elem {
                    return false;
                }
            }
            if one_count != 1 {
                return false;
            }
        }
        true
    }

    pub fn monomial_check(&self) -> bool {
        // Matches the Python reference: every entry is a monomial and sde == 0.
        self.entries.iter().all(|e| e.is_monomial()) && self.sde() == 0
    }

    pub fn hash_value(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }
}

impl PartialEq for Operator {
    fn eq(&self, other: &Self) -> bool {
        self.dim == other.dim
            && self.m == other.m
            && self.n == other.n
            && self.entries == other.entries
    }
}

impl Eq for Operator {}

impl Hash for Operator {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.dim.hash(state);
        self.m.hash(state);
        self.n.hash(state);
        for e in &self.entries {
            e.hash(state);
        }
    }
}
