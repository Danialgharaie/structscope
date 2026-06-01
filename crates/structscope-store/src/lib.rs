use anyhow::Result;
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use structscope_features::FeatureRecord;

mod parquet_out;
pub use parquet_out::write_parquet;

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

    let parquet_path = out_dir.join("features.parquet");
    write_parquet(&parquet_path, records)?;

    let manifest = RunManifest {
        feature_records_path: features_path.display().to_string(),
        parquet_path: Some(parquet_path.display().to_string()),
        notes: vec![
            "Parquet is the primary analytical output; query it with the duckdb-enabled build.".to_string(),
            "JSONL is written alongside as a line-oriented interchange format.".to_string(),
        ],
    };

    let manifest_path = out_dir.join("manifest.json");
    fs::write(&manifest_path, serde_json::to_vec_pretty(&manifest)?)?;
    Ok(manifest)
}

/// Run a SQL query against a feature Parquet file (or directory of them).
/// The data is exposed as a `features` table/view.
pub fn run_query(input: &Path, sql: &str) -> Result<String> {
    #[cfg(feature = "duckdb")]
    {
        query_duckdb(input, sql)
    }

    #[cfg(not(feature = "duckdb"))]
    {
        let _ = (input, sql);
        anyhow::bail!("query support requires building with --features duckdb")
    }
}

#[cfg(feature = "duckdb")]
fn query_duckdb(input: &Path, sql: &str) -> Result<String> {
    use anyhow::Context;
    use duckdb::arrow::util::pretty::pretty_format_batches;
    use duckdb::Connection;

    let target = if input.is_dir() {
        input.join("features.parquet")
    } else {
        input.to_path_buf()
    };
    let target = target.to_str().context("non-UTF8 input path")?;

    let conn = Connection::open_in_memory()?;
    conn.execute_batch(&format!(
        "CREATE VIEW features AS SELECT * FROM read_parquet('{}');",
        target.replace('\'', "''")
    ))?;

    let mut stmt = conn.prepare(sql)?;
    let batches: Vec<_> = stmt.query_arrow([])?.collect();
    Ok(pretty_format_batches(&batches)?.to_string())
}

pub fn normalize_output_path(path: Option<PathBuf>, fallback: &str) -> PathBuf {
    path.unwrap_or_else(|| PathBuf::from(fallback))
}
