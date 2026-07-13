/// Supported cyclotomic dimensions.
/// - 3, 5, 7 are odd-prime qudit dimensions with localization λ = g_p, the
///   quadratic Gauss sum (g_p = √p for p ≡ 1 mod 4, i√p for p ≡ 3 mod 4).
/// - 8 is the qubit case Z[ζ_8, 1/√2].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RingDim {
    D3 = 3,
    D5 = 5,
    D7 = 7,
    D8 = 8,
}

impl RingDim {
    pub fn from_usize(d: usize) -> Option<Self> {
        match d {
            3 => Some(RingDim::D3),
            5 => Some(RingDim::D5),
            7 => Some(RingDim::D7),
            8 => Some(RingDim::D8),
            _ => None,
        }
    }

    pub fn value(self) -> usize {
        self as usize
    }

    /// Dense index into the per-dimension static tables below.
    pub const fn index(self) -> usize {
        match self {
            RingDim::D3 => 0,
            RingDim::D5 => 1,
            RingDim::D7 => 2,
            RingDim::D8 => 3,
        }
    }

    /// Whether this is the qubit (p=8) case, which uses a special reduction branch.
    pub fn is_qubit(self) -> bool {
        matches!(self, RingDim::D8)
    }

    /// Whether complex conjugation flips the sign of the localization:
    /// conj(g_p) = (−1|p)·g_p, so for p ≡ 3 (mod 4) the purely imaginary
    /// λ = i√p conjugates to −λ and the value of α/λ^sde picks up (−1)^sde.
    pub fn conj_flips_sign(self) -> bool {
        matches!(self, RingDim::D3 | RingDim::D7)
    }
}

/// Count how many times each residue is a square modulo p — the coefficient
/// vector of the Gauss sum g_p = Σ ζ^{a²}.
const fn gauss_sequence_const(p: usize) -> [i64; 8] {
    let mut seq = [0i64; 8];
    let mut a = 0;
    while a < p {
        seq[(a * a) % p] += 1;
        a += 1;
    }
    seq
}

/// The `loc_char` matrix used to test divisibility by λ.
///
/// The right-circulant of the gauss sequence (row i = row 0 shifted right by i)
/// represents multiplication by conj(g_p); since g_p·conj(g_p) = p, applying
/// p·circulant and dividing by p² computes exactly α/g_p — valid for every
/// supported p including 8 in principle, though p=8 uses its own branch.
const fn loc_char_const(p: usize) -> [[i64; 8]; 8] {
    let row = gauss_sequence_const(p);
    let mut m = [[0i64; 8]; 8];
    let mut i = 0;
    while i < p {
        let mut j = 0;
        while j < p {
            m[i][j] = (p as i64) * row[(j + p - i) % p];
            j += 1;
        }
        i += 1;
    }
    m
}

/// Compile-time tables, indexed by `RingDim::index()`. Reduction runs on every
/// element construction, so these must never be rebuilt at runtime.
static LOC_CHAR: [[[i64; 8]; 8]; 4] = [
    loc_char_const(3),
    loc_char_const(5),
    loc_char_const(7),
    loc_char_const(8),
];

/// Coefficient vector used to align SDEs during addition: the localization λ
/// itself. For odd primes this is the Gauss sum; for p=8 it is
/// √2 = ζ + ζ⁷ = [0,1,0,0,0,0,0,1] (the literal gauss sequence for n=8 is 4ζ,
/// which is not the localization).
static ALIGNMENT: [[i64; 8]; 4] = [
    gauss_sequence_const(3),
    gauss_sequence_const(5),
    gauss_sequence_const(7),
    [0, 1, 0, 0, 0, 0, 0, 1],
];

pub fn gauss_sequence(dim: RingDim) -> [i64; 8] {
    gauss_sequence_const(dim.value())
}

pub fn loc_char_matrix(dim: RingDim) -> &'static [[i64; 8]; 8] {
    &LOC_CHAR[dim.index()]
}

pub fn sde_alignment_multiplier(dim: RingDim) -> &'static [i64; 8] {
    &ALIGNMENT[dim.index()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_sequence() {
        assert_eq!(gauss_sequence(RingDim::D3)[0..3], [1, 2, 0]);
        assert_eq!(gauss_sequence(RingDim::D5)[0..5], [1, 2, 0, 0, 2]);
        assert_eq!(gauss_sequence(RingDim::D7)[0..7], [1, 2, 2, 0, 2, 0, 0]);
        assert_eq!(gauss_sequence(RingDim::D8), [2, 4, 0, 0, 2, 0, 0, 0]);
    }

    #[test]
    fn test_loc_char_p3_exact() {
        // 3 · right-circulant([1,2,0])
        let m = loc_char_matrix(RingDim::D3);
        assert_eq!(m[0][0..3], [3, 6, 0]);
        assert_eq!(m[1][0..3], [0, 3, 6]);
        assert_eq!(m[2][0..3], [6, 0, 3]);
    }

    #[test]
    fn test_alignment_is_localization() {
        // For odd primes the alignment vector is the Gauss sum itself;
        // for p=8 it must be √2 = ζ + ζ⁷, NOT the gauss sequence.
        assert_eq!(*sde_alignment_multiplier(RingDim::D5), gauss_sequence(RingDim::D5));
        assert_eq!(*sde_alignment_multiplier(RingDim::D8), [0, 1, 0, 0, 0, 0, 0, 1]);
    }

    #[test]
    fn test_conj_sign() {
        assert!(RingDim::D3.conj_flips_sign());
        assert!(!RingDim::D5.conj_flips_sign());
        assert!(RingDim::D7.conj_flips_sign());
        assert!(!RingDim::D8.conj_flips_sign());
    }
}
