use crate::element::CyclotomicElement;
use crate::operator::Operator;
use rustc_hash::FxHashSet;

/// Breadth-first closure of `generators` under multiplication up to `depth`
/// layers. FxHash is deterministic, so iteration order is stable across runs.
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

pub fn torus(subgroup: &[Operator], null_elem: &CyclotomicElement) -> Vec<Operator> {
    subgroup
        .iter()
        .filter(|op| op.is_diag(null_elem))
        .cloned()
        .collect()
}

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

/// Right (or left) coset representatives of H in G.
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
        reps.push(coset.into_iter().next().unwrap_or(elem));
    }

    reps
}

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

pub fn multiply_many(operators: &[Operator]) -> Option<Operator> {
    let mut iter = operators.iter();
    let mut acc = iter.next()?.clone();
    for op in iter {
        acc = acc.matmul(op);
    }
    Some(acc)
}

/// Prefix products of a left-multiplication walk driven by `indices`:
/// returns [g[i0], g[i1]*g[i0], ..., g[ik]*...*g[i0]]. Batches random-walk
/// sampling into one FFI call — the caller picks the random indices cheaply
/// in Python, the products stay in Rust.
pub fn multiply_selected_chain(generators: &[Operator], indices: &[usize]) -> Option<Vec<Operator>> {
    let mut iter = indices.iter();
    let mut acc = generators.get(*iter.next()?)?.clone();
    let mut out = Vec::with_capacity(indices.len());
    out.push(acc.clone());
    for &i in iter {
        acc = generators.get(i)?.matmul(&acc);
        out.push(acc.clone());
    }
    Some(out)
}

pub fn tensor_power(op: &Operator, power: usize) -> Option<Operator> {
    if power == 0 {
        return None;
    }
    let mut acc = op.clone();
    for _ in 1..power {
        acc = acc.tensor(op);
    }
    Some(acc)
}
