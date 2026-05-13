use anyhow::{bail, Result};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use structscope_features::FeatureRecord;

#[derive(Debug, Clone, Serialize)]
pub struct RunManifest {
    pub feature_records_path: String,
    pub parquet_path: Option<String>,
    pub notes: Vec<String>,
}

pub fn write_feature_records(out_dir: &Path, records: &[FeatureRecord]) -> Result<RunManifest> {
    fs::create_dir_all(out_dir)?;
    let features_path = out_dir.join("features.jsonl");
    let mut file = File::create(&features_path)?;
    for record in records {
        writeln!(file, "{}", serde_json::to_string(record)?)?;
    }

    let manifest = RunManifest {
        feature_records_path: features_path.display().to_string(),
        parquet_path: None,
        notes: vec![
            "Parquet output is not implemented in this bootstrap slice.".to_string(),
            "JSONL output is written as the current analytical interchange format.".to_string(),
        ],
    };

    let manifest_path = out_dir.join("manifest.json");
    fs::write(&manifest_path, serde_json::to_vec_pretty(&manifest)?)?;
    Ok(manifest)
}

pub fn run_query(_input: &Path, _sql: &str) -> Result<String> {
    #[cfg(feature = "duckdb")]
    {
        bail!("DuckDB support is declared but not implemented in this bootstrap slice")
    }

    #[cfg(not(feature = "duckdb"))]
    {
        bail!("query support requires a future DuckDB integration build")
    }
}

pub fn normalize_output_path(path: Option<PathBuf>, fallback: &str) -> PathBuf {
    path.unwrap_or_else(|| PathBuf::from(fallback))
}
