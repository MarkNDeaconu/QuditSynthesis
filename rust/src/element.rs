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

/// Raw cyclic convolution of two elements' coefficient vectors into an i128
/// accumulator (no reduction, no sde bookkeeping).
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

/// Cyclic convolution of an i128 accumulator with an i64 multiplier vector.
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
    /// Create a new element and reduce it to canonical form.
    /// The canonical zero always has sde 0.
    pub fn new(dim: RingDim, coeffs: [i64; 8], sde: i32) -> Self {
        Self::from_i128(dim, widen(&coeffs), sde)
    }

    /// Canonicalize an i128 coefficient vector: divide out the localization
    /// while possible, subtract the mode, normalize zero, then narrow to i64.
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

    /// Add two elements, aligning denominators (multiplying the lower-sde side
    /// by λ per step) in i128, then reduce once.
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

    /// Multiply by another cyclotomic element and reduce.
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.dim, other.dim);
        Self::from_i128(self.dim, conv_pair_i128(self, other), self.sde + other.sde)
    }

    /// Multiply by a scalar integer.
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

    /// Evaluate the complex value of the element, given the localization as (re, im).
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
        // Divide by localization^sde using polar form.
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

/// Divide out λ while possible (divisibility test via the integer loc_char
/// circulant, congruence mod p²), then subtract the mode — all in i128.
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
    // the deterministic rule shared with the Python code). O(p²) scan, no alloc.
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

    // Divide by (√2)^2 = 2 while possible.
    while [a, b, c, d].iter().all(|x| x.rem_euclid(2) == 0) && (a != 0 || b != 0 || c != 0 || d != 0)
    {
        a /= 2;
        b /= 2;
        c /= 2;
        d /= 2;
        *sde -= 2;
    }

    // Divide by √2 if the residue pattern matches.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::RingDim;

    fn e3(coeffs: [i64; 3], sde: i32) -> CyclotomicElement {
        let mut c = [0i64; 8];
        c[..3].copy_from_slice(&coeffs);
        CyclotomicElement::new(RingDim::D3, c, sde)
    }

    #[test]
    fn test_add_same_sde() {
        // 1 + ζ = -ζ²
        let sum = e3([1, 0, 0], 0).add(&e3([0, 1, 0], 0));
        assert_eq!(sum.coeffs[0..3], [0, 0, -1]);
    }

    #[test]
    fn test_mul() {
        // ζ * ζ = ζ²
        let z = e3([0, 1, 0], 0);
        assert_eq!(z.mul(&z).coeffs[0..3], [0, 0, 1]);
    }

    #[test]
    fn test_reduce_sde() {
        // g_3 = 1 + 2ζ = i√3, so [1,2,0] should reduce to [1,0,0] with sde=-1.
        let e = e3([1, 2, 0], 0);
        assert_eq!(e.coeffs[0..3], [1, 0, 0]);
        assert_eq!(e.sde, -1);
    }

    #[test]
    fn test_conj_sign_odd_sde() {
        // e = 1/(i√3): conj must be 1/(-i√3) = -1/(i√3), i.e. numerator -1.
        let e = e3([1, 0, 0], 1);
        let c = e.conj();
        assert_eq!(c.coeffs[0..3], [-1, 0, 0]);
        assert_eq!(c.sde, 1);
        // Involution: conj(conj(e)) == e.
        assert_eq!(c.conj(), e);
        // Even sde: no sign.
        let f = e3([0, 1, 0], 2);
        assert_eq!(f.conj().coeffs[0..3], [0, 0, 1]);
        // Real localization (p=5): never flips.
        let mut c5 = [0i64; 8];
        c5[0] = 1;
        let g = CyclotomicElement::new(RingDim::D5, c5, 1);
        assert_eq!(g.conj().coeffs[0..5], [1, 0, 0, 0, 0]);
    }

    #[test]
    fn test_norm_positive_via_conj() {
        // e·conj(e) for e = 1/(i√3) must be +1/3, i.e. [−1,0,0] at sde 2
        // (1/3 = −1/λ² since λ² = −3).
        let e = e3([1, 0, 0], 1);
        let n = e.mul(&e.conj());
        assert_eq!(n.coeffs[0..3], [-1, 0, 0]);
        assert_eq!(n.sde, 2);
        let loc = (0.0, 3f64.sqrt());
        let (re, im) = n.comp(loc);
        assert!((re - 1.0 / 3.0).abs() < 1e-12 && im.abs() < 1e-12);
    }

    #[test]
    fn test_disguised_zero_canonicalizes() {
        // [1,1,1] = 1+ζ+ζ² = 0 must canonicalize to sde 0 regardless of input sde.
        let z = e3([1, 1, 1], 5);
        assert!(z.is_zero());
        assert_eq!(z.sde, 0);
        assert_eq!(z, e3([0, 0, 0], 0));
    }

    #[test]
    fn test_mode_tiebreak_largest_value() {
        // All counts equal → subtract the largest value (deterministic rule).
        // [5,7,1] has coefficient sum ≢ 0 (mod 3) so it is not divisible by g_3;
        // counts all 1, mode = 7 → [-2,0,-6].
        let e = e3([5, 7, 1], 0);
        assert_eq!(e.coeffs[0..3], [-2, 0, -6]);
    }

    #[test]
    fn test_scale_localization() {
        // p=5: (1/√5)·1 → sde −(−1)... i.e. [1,...] at sde+1.
        let mut c5 = [0i64; 8];
        c5[0] = 1;
        let e = CyclotomicElement::new(RingDim::D5, c5, 0);
        let s = e.scale_localization(1, -1).unwrap();
        assert_eq!((s.coeffs[0], s.sde), (1, 1));
        // p=3: 1/3 = |λ|^{-2} with sign flip: −coeffs at sde+2.
        let e = e3([0, 1, 0], 0);
        let s = e.scale_localization(1, -2).unwrap();
        assert_eq!(s.coeffs[0..3], [0, -1, 0]);
        assert_eq!(s.sde, 2);
        // p=3: odd powers of √3 are not ring elements.
        assert!(e.scale_localization(1, -1).is_err());
    }

    #[test]
    fn test_add_large_sde_gap_exact() {
        // A big sde gap forces long alignment chains whose intermediates exceed
        // i64; the i128 pipeline must still produce the exact canonical result.
        let a = e3([1, 0, 0], 30);
        let b = e3([0, 1, 0], 0);
        let s = a.add(&b);
        // Subtracting a back must recover b exactly.
        let neg_a = a.mul_scalar(-1);
        assert_eq!(s.add(&neg_a), b);
    }
}
