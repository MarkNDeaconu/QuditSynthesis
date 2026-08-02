use crate::ring::{loc_char_matrix, sde_alignment_multiplier, RingDim};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// A cyclotomic element in Z[ζ_p, 1/λ] represented as a coefficient vector
/// of length `dim.value()` plus a smallest-denominator exponent (sde),
/// with value α/λ^sde.
///
/// Canonical coefficients are stored in i64, but all arithmetic accumulates
/// and reduces in i128: intermediate sums (raw matmul accumulators, sde
/// alignment products) can legitimately exceed i64 even when the canonical
/// result fits. Only the final narrowing can overflow, and it panics with a
/// clear message (PyO3 converts that into a catchable PanicException).
#[derive(Clone, Debug)]
pub struct CyclotomicElement {
    pub dim: RingDim,
    pub coeffs: [i64; 8],
    pub sde: i32,
}

pub(crate) fn conv_pair_i128(a: &CyclotomicElement, b: &CyclotomicElement) -> [i128; 8] {
    let p = a.dimension();
    let mut out = [0i128; 8];
    for i in 0..p {
        if a.coeffs[i] == 0 {
            continue;
        }
        for j in 0..p {
            if b.coeffs[j] == 0 {
                continue;
            }
            out[(i + j) % p] += a.coeffs[i] as i128 * b.coeffs[j] as i128;
        }
    }
    out
}

pub(crate) fn conv_i128(a: &[i128; 8], b: &[i64; 8], p: usize) -> [i128; 8] {
    let mut out = [0i128; 8];
    for i in 0..p {
        if a[i] == 0 {
            continue;
        }
        for j in 0..p {
            if b[j] == 0 {
                continue;
            }
            out[(i + j) % p] += a[i] * b[j] as i128;
        }
    }
    out
}

fn widen(coeffs: &[i64; 8]) -> [i128; 8] {
    let mut out = [0i128; 8];
    for i in 0..8 {
        out[i] = coeffs[i] as i128;
    }
    out
}

impl CyclotomicElement {
    /// Reduce to canonical form on construction; the canonical zero always has sde 0.
    pub fn new(dim: RingDim, coeffs: [i64; 8], sde: i32) -> Self {
        Self::from_i128(dim, widen(&coeffs), sde)
    }

    pub(crate) fn from_i128(dim: RingDim, mut coeffs: [i128; 8], mut sde: i32) -> Self {
        let p = dim.value();
        if coeffs[0..p].iter().any(|&x| x != 0) {
            if dim.is_qubit() {
                reduce_qubit_i128(&mut coeffs, &mut sde);
            } else {
                reduce_prime_i128(dim, &mut coeffs, &mut sde);
            }
        }
        let mut out = [0i64; 8];
        let mut is_zero = true;
        for i in 0..p {
            out[i] = coeffs[i]
                .try_into()
                .expect("canonical coefficient overflow: value exceeds i64");
            if out[i] != 0 {
                is_zero = false;
            }
        }
        Self {
            dim,
            coeffs: out,
            sde: if is_zero { 0 } else { sde },
        }
    }

    pub fn dimension(&self) -> usize {
        self.dim.value()
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs[0..self.dimension()].iter().all(|x| *x == 0)
    }

    pub fn is_monomial(&self) -> bool {
        self.coeffs[0..self.dimension()]
            .iter()
            .filter(|x| **x != 0)
            .count()
            <= 1
    }

    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        let p = self.dimension();
        let denom = sde_alignment_multiplier(self.dim);
        let target = self.sde.max(other.sde);
        let mut acc = widen(&self.coeffs);
        for _ in 0..(target - self.sde) {
            acc = conv_i128(&acc, denom, p);
        }
        let mut rhs = widen(&other.coeffs);
        for _ in 0..(target - other.sde) {
            rhs = conv_i128(&rhs, denom, p);
        }
        for i in 0..p {
            acc[i] += rhs[i];
        }
        Self::from_i128(self.dim, acc, target)
    }

    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        Self::from_i128(self.dim, conv_pair_i128(self, other), self.sde + other.sde)
    }

    pub fn mul_scalar(&self, scalar: i64) -> Self {
        let p = self.dimension();
        let mut acc = [0i128; 8];
        for i in 0..p {
            acc[i] = self.coeffs[i] as i128 * scalar as i128;
        }
        Self::from_i128(self.dim, acc, self.sde)
    }

    /// Multiply by the exact scalar sign·|λ|^k (e.g. 1/√p is sign=1, k=−1),
    /// absorbing the power into the sde. For p ≡ 3 (mod 4), |λ|^k = p^{k/2}
    /// requires k even (√p itself is not in the ring) and equals (−1)^{k/2}·λ^k,
    /// contributing a sign; for real λ (p = 5, 8) the power absorbs directly.
    pub fn scale_localization(&self, sign: i64, k: i32) -> Result<Self, &'static str> {
        let mut sign = sign;
        if self.dim.conj_flips_sign() {
            if k.rem_euclid(2) == 1 {
                return Err("odd powers of √p are not ring elements for p ≡ 3 (mod 4)");
            }
            if k.div_euclid(2).rem_euclid(2) == 1 {
                sign = -sign;
            }
        }
        let mut coeffs = self.coeffs;
        if sign < 0 {
            for c in coeffs[0..self.dimension()].iter_mut() {
                *c = -*c;
            }
        }
        Ok(Self::new(self.dim, coeffs, self.sde - k))
    }

    /// Complex conjugation: ζ ↦ ζ^{-1} and λ ↦ conj(λ). For p ≡ 3 (mod 4) the
    /// localization is purely imaginary (conj(g_p) = −g_p), so the value picks
    /// up a factor (−1)^sde which must land in the numerator.
    pub fn conj(&self) -> Self {
        let p = self.dimension();
        let mut coeffs = [0i64; 8];
        coeffs[0] = self.coeffs[0];
        for i in 1..p {
            coeffs[i] = self.coeffs[p - i];
        }
        if self.dim.conj_flips_sign() && self.sde.rem_euclid(2) == 1 {
            for c in coeffs[0..p].iter_mut() {
                *c = -*c;
            }
        }
        Self::new(self.dim, coeffs, self.sde)
    }

    pub fn comp(&self, localization: (f64, f64)) -> (f64, f64) {
        let p = self.dimension();
        let theta = 2.0 * std::f64::consts::PI / (p as f64);
        let mut re = 0.0;
        let mut im = 0.0;
        for i in 0..p {
            let (zr, zi) = ((theta * i as f64).cos(), (theta * i as f64).sin());
            re += self.coeffs[i] as f64 * zr;
            im += self.coeffs[i] as f64 * zi;
        }
        let (lr, li) = localization;
        let r = (lr * lr + li * li).sqrt();
        let arg = li.atan2(lr);
        let sde = self.sde as f64;
        let denom_r = r.powf(sde);
        let denom_arg = arg * sde;
        let denom_re = denom_r * denom_arg.cos();
        let denom_im = denom_r * denom_arg.sin();
        let denom_norm = denom_re * denom_re + denom_im * denom_im;
        if denom_norm == 0.0 {
            return (f64::INFINITY, f64::INFINITY);
        }
        let nr = re * denom_re + im * denom_im;
        let ni = im * denom_re - re * denom_im;
        (nr / denom_norm, ni / denom_norm)
    }

    pub fn hash_value(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }
}

/// Divide out λ while possible (loc_char circulant divisibility test, congruence
/// mod p²), then subtract the mode — all in i128.
fn reduce_prime_i128(dim: RingDim, coeffs: &mut [i128; 8], sde: &mut i32) {
    let p = dim.value();
    let p2 = (p * p) as i128;
    let loc = loc_char_matrix(dim);

    loop {
        // If all coefficients are equal, the element is k*(1+ζ+...+ζ^{p-1}) = 0.
        if coeffs[0..p].iter().all(|&x| x == coeffs[0]) {
            break;
        }
        let mut nc = [0i128; 8];
        for i in 0..p {
            for j in 0..p {
                nc[i] += loc[i][j] as i128 * coeffs[j];
            }
        }
        let r = nc[0].rem_euclid(p2);
        if nc[0..p].iter().all(|&x| x.rem_euclid(p2) == r) {
            for i in 0..p {
                coeffs[i] = (nc[i] - r) / p2;
            }
            *sde -= 1;
        } else {
            break;
        }
    }

    // Subtract the most frequent coefficient (ties break by the largest value —
    // the deterministic rule shared with the Python code).
    let mut best = coeffs[0];
    let mut best_count = 0usize;
    for i in 0..p {
        let v = coeffs[i];
        let mut count = 0usize;
        for j in 0..p {
            if coeffs[j] == v {
                count += 1;
            }
        }
        if count > best_count || (count == best_count && v > best) {
            best = v;
            best_count = count;
        }
    }
    if best != 0 {
        for c in coeffs[0..p].iter_mut() {
            *c -= best;
        }
    }
}

/// The p=8 qubit reduction: fold to 4 coefficients via ζ⁴ = −1, divide by 2
/// while possible, then by √2 when the residue pattern allows — all in i128.
fn reduce_qubit_i128(coeffs: &mut [i128; 8], sde: &mut i32) {
    let mut a = coeffs[0] - coeffs[4];
    let mut b = coeffs[1] - coeffs[5];
    let mut c = coeffs[2] - coeffs[6];
    let mut d = coeffs[3] - coeffs[7];

    while [a, b, c, d].iter().all(|x| x.rem_euclid(2) == 0) && (a != 0 || b != 0 || c != 0 || d != 0)
    {
        a /= 2;
        b /= 2;
        c /= 2;
        d /= 2;
        *sde -= 2;
    }

    let r = [a.rem_euclid(2), b.rem_euclid(2), c.rem_euclid(2), d.rem_euclid(2)];
    if r == [1, 0, 1, 0] || r == [0, 1, 0, 1] || r == [1, 1, 1, 1] {
        let (a2, b2, c2, d2) = ((b - d) / 2, (c + a) / 2, (b + d) / 2, (c - a) / 2);
        a = a2;
        b = b2;
        c = c2;
        d = d2;
        *sde -= 1;
    }

    *coeffs = [a, b, c, d, 0, 0, 0, 0];
}

impl PartialEq for CyclotomicElement {
    fn eq(&self, other: &Self) -> bool {
        self.dim == other.dim
            && self.sde == other.sde
            && self.coeffs[0..self.dimension()] == other.coeffs[0..other.dimension()]
    }
}

impl Eq for CyclotomicElement {}

impl Hash for CyclotomicElement {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.dim.hash(state);
        self.sde.hash(state);
        self.coeffs[0..self.dimension()].hash(state);
    }
}
