//! Global sequence alignment (Needleman-Wunsch) for establishing residue
//! correspondence between two structures. Raw primitive; emits matched index pairs.

/// One-letter code for a standard amino acid three-letter name, else `b'X'`.
pub fn three_to_one(residue_name: &str) -> u8 {
    match residue_name.trim().to_ascii_uppercase().as_str() {
        "ALA" => b'A',
        "ARG" => b'R',
        "ASN" => b'N',
        "ASP" => b'D',
        "CYS" => b'C',
        "GLN" => b'Q',
        "GLU" => b'E',
        "GLY" => b'G',
        "HIS" => b'H',
        "ILE" => b'I',
        "LEU" => b'L',
        "LYS" => b'K',
        "MET" => b'M',
        "PHE" => b'F',
        "PRO" => b'P',
        "SER" => b'S',
        "THR" => b'T',
        "TRP" => b'W',
        "TYR" => b'Y',
        "VAL" => b'V',
        _ => b'X',
    }
}

/// Global-align `a` and `b` (identity scoring: match +1, mismatch -1, gap -1) and
/// return the aligned column index pairs `(i, j)` where neither side is a gap.
pub fn needleman_wunsch(a: &[u8], b: &[u8]) -> Vec<(usize, usize)> {
    let (n, m) = (a.len(), b.len());
    const GAP: i32 = -1;
    let mut score = vec![vec![0i32; m + 1]; n + 1];
    for i in 0..=n {
        score[i][0] = GAP * i as i32;
    }
    for j in 0..=m {
        score[0][j] = GAP * j as i32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let s = if a[i - 1] == b[j - 1] { 1 } else { -1 };
            score[i][j] = (score[i - 1][j - 1] + s).max(score[i - 1][j] + GAP).max(score[i][j - 1] + GAP);
        }
    }

    let mut pairs = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        let s = if a[i - 1] == b[j - 1] { 1 } else { -1 };
        if score[i][j] == score[i - 1][j - 1] + s {
            pairs.push((i - 1, j - 1));
            i -= 1;
            j -= 1;
        } else if score[i][j] == score[i - 1][j] + GAP {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    pairs.reverse();
    pairs
}

/// Local-align `a` and `b` (Smith-Waterman: match +2, mismatch -1, gap -1) and
/// return the matched index pairs `(i, j)` of the single best local region.
/// Useful for partial/domain overlaps where global alignment is inappropriate.
pub fn smith_waterman(a: &[u8], b: &[u8]) -> Vec<(usize, usize)> {
    let (n, m) = (a.len(), b.len());
    const GAP: i32 = -1;
    let mut score = vec![vec![0i32; m + 1]; n + 1];
    let (mut best, mut best_i, mut best_j) = (0i32, 0usize, 0usize);
    for i in 1..=n {
        for j in 1..=m {
            let s = if a[i - 1] == b[j - 1] { 2 } else { -1 };
            let v = (score[i - 1][j - 1] + s).max(score[i - 1][j] + GAP).max(score[i][j - 1] + GAP).max(0);
            score[i][j] = v;
            if v > best {
                best = v;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Trace back from the maximum cell until a zero score is reached.
    let mut pairs = Vec::new();
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 && score[i][j] > 0 {
        let s = if a[i - 1] == b[j - 1] { 2 } else { -1 };
        if score[i][j] == score[i - 1][j - 1] + s {
            pairs.push((i - 1, j - 1));
            i -= 1;
            j -= 1;
        } else if score[i][j] == score[i - 1][j] + GAP {
            i -= 1;
        } else {
            j -= 1;
        }
    }
    pairs.reverse();
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_align_fully() {
        let a = b"ACDEFG";
        let pairs = needleman_wunsch(a, a);
        assert_eq!(pairs.len(), 6);
        assert!(pairs.iter().enumerate().all(|(k, &(i, j))| i == k && j == k));
    }

    #[test]
    fn handles_internal_gap() {
        // b is a with one residue deleted: 5 columns should still match.
        let pairs = needleman_wunsch(b"ACDEFG", b"ACEFG");
        assert_eq!(pairs.len(), 5);
        // First two positions align 1:1.
        assert_eq!(pairs[0], (0, 0));
        assert_eq!(pairs[1], (1, 1));
    }

    #[test]
    fn three_to_one_maps_standard_and_unknown() {
        assert_eq!(three_to_one("ALA"), b'A');
        assert_eq!(three_to_one("trp"), b'W');
        assert_eq!(three_to_one("HOH"), b'X');
    }

    #[test]
    fn smith_waterman_finds_local_region() {
        // Shared core "DEFGH" flanked by dissimilar sequence on both sides.
        let pairs = smith_waterman(b"WWWDEFGHWW", b"YYDEFGHYYY");
        assert_eq!(pairs.len(), 5);
        // The matched columns map the shared region 1:1.
        let a = b"WWWDEFGHWW";
        let b = b"YYDEFGHYYY";
        assert!(pairs.iter().all(|&(i, j)| a[i] == b[j]));
        assert_eq!(pairs[0], (3, 2));
    }
}
