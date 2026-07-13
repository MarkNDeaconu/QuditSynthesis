use crate::element::CyclotomicElement;
use crate::operator::Operator;
use rustc_hash::FxHashSet;

/// Breadth-first closure of `generators` under multiplication up to `depth` layers.
/// Frontier-based: each layer only multiplies generators against the previous
/// layer's new elements, with early exit once the closure stabilizes.
/// (FxHash is deterministic, so iteration order is stable across runs.)
pub fn subgroup_bfs(generators: &[Operator], depth: usize) -> Vec<Operator> {
    if generators.is_empty() {
        return Vec::new();
    }
    let mut orbit: FxHashSet<Operator> = generators.iter().cloned().collect();
    let mut current: Vec<Operator> = orbit.iter().cloned().collect();

    for _ in 0..depth {
        let mut next = Vec::new();
        for g in generators {
            for o in &current {
                let prod = g.matmul(o);
                if !orbit.contains(&prod) {
                    next.push(prod.clone());
                    orbit.insert(prod);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        current = next;
    }

    orbit.into_iter().collect()
}

/// Filter `subgroup` to diagonal operators.
pub fn torus(subgroup: &[Operator], null_elem: &CyclotomicElement) -> Vec<Operator> {
    subgroup
        .iter()
        .filter(|op| op.is_diag(null_elem))
        .cloned()
        .collect()
}

/// Filter `subgroup` to permutation operators (exactly one `one_elem` per row,
/// all other entries equal `null_elem`).
pub fn permutation_subgroup(
    subgroup: &[Operator],
    one_elem: &CyclotomicElement,
    null_elem: &CyclotomicElement,
) -> Vec<Operator> {
    subgroup
        .iter()
        .filter(|op| op.is_permutation(one_elem, null_elem))
        .cloned()
        .collect()
}

/// Compute a set of right (or left) coset representatives of H in G.
pub fn quotient(g: &[Operator], h: &[Operator], right: bool) -> Vec<Operator> {
    let mut group: FxHashSet<Operator> = g.iter().cloned().collect();
    let mut reps = Vec::new();

    while let Some(elem) = group.iter().next().cloned() {
        group.remove(&elem);
        let coset: FxHashSet<Operator> = if right {
            h.iter().map(|hh| elem.matmul(hh)).collect()
        } else {
            h.iter().map(|hh| hh.matmul(&elem)).collect()
        };
        for c in &coset {
            group.remove(c);
        }
        // Pop a representative from the coset (deterministic with FxHash).
        reps.push(coset.into_iter().next().unwrap_or(elem));
    }

    reps
}

/// Left-multiply `op` by each element of `dropping_set` and return the first
/// product whose total SDE is strictly smaller, together with the gate string.
pub fn synth_search(op: &Operator, dropping_set: &[Operator]) -> Option<(Operator, String)> {
    let base = op.sde_sum();
    for gate in dropping_set {
        let candidate = gate.matmul(op);
        if candidate.sde_sum() < base {
            return Some((candidate, gate.gate_string.clone()));
        }
    }
    None
}

/// Repeatedly apply `synth_search` until the minimum entry SDE is at most
/// `target_sde`. Returns the concatenated gate string, or `None` if the search
/// stalls (the PyO3 layer converts that into a RuntimeError, matching the
/// Python reference).
pub fn synthesize(op: &Operator, dropping_set: &[Operator], target_sde: i32) -> Option<String> {
    let mut mat = op.clone();
    // Each step's gate string is PREPENDED to the result; collect and join once.
    let mut labels: Vec<String> = Vec::new();

    loop {
        let min_sde = mat
            .entries
            .iter()
            .map(|e| e.sde)
            .min()
            .unwrap_or(i32::MIN);
        if min_sde <= target_sde {
            let mut out = String::with_capacity(labels.iter().map(|s| s.len()).sum());
            for s in labels.iter().rev() {
                out.push_str(s);
            }
            return Some(out);
        }

        let (new_mat, string) = synth_search(&mat, dropping_set)?;
        mat = new_mat;
        labels.push(string);
    }
}

/// Multiply a list of operators left-to-right, keeping the whole product in Rust.
pub fn multiply_many(operators: &[Operator]) -> Option<Operator> {
    let mut iter = operators.iter();
    let mut acc = iter.next()?.clone();
    for op in iter {
        acc = acc.matmul(op);
    }
    Some(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operator::Operator;
    use crate::ring::RingDim;

    #[test]
    fn test_subgroup_bfs_identity() {
        let id = Operator::identity(RingDim::D3, 2);
        let gens = vec![id.clone()];
        let orbit = subgroup_bfs(&gens, 3);
        assert_eq!(orbit.len(), 1);
        assert_eq!(orbit[0], id);
    }

    #[test]
    fn test_quotient_self() {
        let id = Operator::identity(RingDim::D3, 2);
        let g = vec![id.clone()];
        let reps = quotient(&g, &g, true);
        assert_eq!(reps.len(), 1);
    }

    #[test]
    fn test_multiply_many_identity() {
        let id = Operator::identity(RingDim::D3, 2);
        let result = multiply_many(&[id.clone(), id.clone(), id.clone()]).unwrap();
        assert_eq!(result, id);
    }

    #[test]
    fn test_synthesize_prepend_order() {
        // The concatenated string must be s_last + ... + s_first, matching the
        // reference's `final_string = string + final_string` per step.
        let d = RingDim::D3;
        let mut c = [0i64; 8];
        c[0] = 1;
        let one = crate::element::CyclotomicElement::new(d, c, 0);
        let mut c2 = [0i64; 8];
        c2[..3].copy_from_slice(&[1, 2, 0]); // g_3: multiplying by it lowers sde by 1... (sde -1 canonical)
        let g_elem = crate::element::CyclotomicElement::new(d, c2, 0);
        // target = diag(g^4-ish high sde element): build diag with sde 3 via [1,0,0] sde 3
        let mut c3 = [0i64; 8];
        c3[0] = 1;
        let high = crate::element::CyclotomicElement::new(d, c3, 3);
        let zero = crate::element::CyclotomicElement::new(d, [0i64; 8], 0);
        let target = Operator::new(d, 1, 1, vec![high], String::new());
        let mut drop_gate = Operator::new(d, 1, 1, vec![g_elem], String::new());
        drop_gate.gate_string = "G".to_string();
        let _ = (one, zero);
        let s = synthesize(&target, &[drop_gate], 1).unwrap();
        assert_eq!(s, "GG"); // two drops: sde 3 -> 2 -> 1
    }
}
