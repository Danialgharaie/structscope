use crate::Residue;

const DEFAULT_DENYLIST: &[&str] = &[
    "HOH", "WAT", "DOD", "NA", "CL", "MG", "ZN", "CA", "K", "FE", "MN", "CU", "CO", "CD", "NI",
    "SR", "BA", "CS", "RB", "LI", "IOD",
];

#[derive(Debug, Clone, Default)]
pub struct LigandFilter {
    include_only: Option<Vec<String>>,
    extra_exclude: Vec<String>,
}

impl LigandFilter {
    pub fn include_only(names: &[&str]) -> Self {
        Self {
            include_only: Some(names.iter().map(|s| s.to_ascii_uppercase()).collect()),
            extra_exclude: vec![],
        }
    }

    pub fn with_extra_exclude(mut self, names: &[&str]) -> Self {
        self.extra_exclude
            .extend(names.iter().map(|s| s.to_ascii_uppercase()));
        self
    }

    pub fn is_ligand(&self, residue: &Residue) -> bool {
        if !residue.is_hetero {
            return false;
        }
        let name = residue.name.trim().to_ascii_uppercase();
        if let Some(allow) = &self.include_only {
            return allow.iter().any(|n| n == &name);
        }
        if DEFAULT_DENYLIST.iter().any(|d| d.eq_ignore_ascii_case(&name)) {
            return false;
        }
        !self.extra_exclude.iter().any(|n| n == &name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Residue;

    fn het_res(name: &str) -> Residue {
        Residue {
            id: format!("A:{name}:1"),
            name: name.to_string(),
            seq_number: 1,
            insertion_code: None,
            atoms: vec![],
            is_hetero: true,
        }
    }

    fn std_res(name: &str) -> Residue {
        let mut r = het_res(name);
        r.is_hetero = false;
        r
    }

    #[test]
    fn default_excludes_water_and_sodium() {
        let f = LigandFilter::default();
        assert!(!f.is_ligand(&het_res("HOH")));
        assert!(!f.is_ligand(&het_res("NA")));
        assert!(f.is_ligand(&het_res("HEM")));
    }

    #[test]
    fn exclude_is_additive() {
        let f = LigandFilter::default().with_extra_exclude(&["HEM"]);
        assert!(!f.is_ligand(&het_res("HEM")));
        assert!(f.is_ligand(&het_res("NAG")));
    }

    #[test]
    fn include_only_overrides_defaults() {
        let f = LigandFilter::include_only(&["NA"]);
        assert!(f.is_ligand(&het_res("NA")));
        assert!(!f.is_ligand(&het_res("HEM")));
    }

    #[test]
    fn protein_residues_never_ligands() {
        let f = LigandFilter::default();
        assert!(!f.is_ligand(&std_res("ALA")));
    }
}
