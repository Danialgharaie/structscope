use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use serde_json::json;
use std::fs;
use std::path::{Path, PathBuf};
use structscope_agent::guard_available;
use structscope_core::{parse_file, ParseOptions, Structure};
use structscope_events::Event;
use structscope_features::compute_features;
use structscope_graphs::{build_residue_graph, export_graphml};
use structscope_provenance::{inspect_sqlite, ProvenanceConfig, ProvenanceRecorder};
use structscope_store::{normalize_output_path, run_query, write_feature_records};
use walkdir::WalkDir;

#[derive(Parser)]
#[command(name = "structscope", version, about = "Structural bioinformatics toolkit bootstrap")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Parse {
        input: PathBuf,
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        format: OutputFormat,
    },
    Featurize {
        input: PathBuf,
        #[arg(long, default_value = "all")]
        features: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        provenance: bool,
        #[arg(long)]
        sqlite: Option<PathBuf>,
        #[arg(long)]
        jsonl: Option<PathBuf>,
        #[arg(long)]
        guard: bool,
    },
    Graph {
        input: PathBuf,
        #[arg(long, default_value = "residue")]
        graph_type: String,
        #[arg(long, default_value = "graphml")]
        format: String,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    Query {
        input: PathBuf,
        #[arg(long)]
        sql: String,
    },
    Provenance {
        sqlite: PathBuf,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OutputFormat {
    Text,
    Json,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Parse { input, format } => cmd_parse(&input, format),
        Commands::Featurize {
            input,
            features: _,
            out,
            provenance,
            sqlite,
            jsonl,
            guard,
        } => cmd_featurize(&input, &out, provenance, sqlite, jsonl, guard),
        Commands::Graph {
            input,
            graph_type,
            format,
            out,
        } => cmd_graph(&input, &graph_type, &format, out),
        Commands::Query { input, sql } => {
            println!("{}", run_query(&input, &sql)?);
            Ok(())
        }
        Commands::Provenance { sqlite } => cmd_provenance(&sqlite),
    }
}

fn cmd_parse(input: &Path, format: OutputFormat) -> Result<()> {
    let structure = parse_file(input, ParseOptions::default())?;
    match format {
        OutputFormat::Text => {
            let summary = structure.summary();
            println!(
                "structure_id={}; chains={}; residues={}; atoms={}; heteroatoms={}; ligands={}",
                summary.structure_id,
                summary.chain_count,
                summary.residue_count,
                summary.atom_count,
                summary.heteroatom_count,
                summary.ligand_count
            );
        }
        OutputFormat::Json => {
            println!("{}", serde_json::to_string_pretty(&structure.summary())?);
        }
    }
    Ok(())
}

fn cmd_featurize(
    input: &Path,
    out: &Path,
    provenance: bool,
    sqlite: Option<PathBuf>,
    jsonl: Option<PathBuf>,
    guard: bool,
) -> Result<()> {
    let inputs = collect_inputs(input)?;
    fs::create_dir_all(out)?;
    let mut recorder = if provenance {
        Some(ProvenanceRecorder::open(
            &ProvenanceConfig {
                sqlite_path: Some(normalize_output_path(
                    sqlite,
                    &out.join("run.sqlite").display().to_string(),
                )),
                jsonl_path: Some(normalize_output_path(
                    jsonl,
                    &out.join("events.jsonl").display().to_string(),
                )),
            },
            "featurize",
        )?)
    } else {
        None
    };

    if guard && !guard_available() {
        eprintln!("guard requested, but the optional eBPF agent is not implemented in this bootstrap slice");
    }

    let mut records = Vec::new();
    for path in inputs {
        if let Some(recorder) = recorder.as_mut() {
            recorder.record(Event::new(
                "parse_start",
                None,
                json!({ "path": path.display().to_string() }),
            ))?;
        }

        match parse_file(&path, ParseOptions::default()) {
            Ok(structure) => {
                let record = compute_features(&structure);
                if let Some(recorder) = recorder.as_mut() {
                    recorder.record(Event::new(
                        "feature_complete",
                        Some(structure.id.clone()),
                        json!({ "path": path.display().to_string() }),
                    ))?;
                }
                records.push(record);
            }
            Err(err) => {
                if let Some(recorder) = recorder.as_mut() {
                    recorder.record(Event::new(
                        "structure_failed",
                        None,
                        json!({ "path": path.display().to_string(), "error": err.to_string() }),
                    ))?;
                }
            }
        }
    }

    let manifest = write_feature_records(out, &records)?;
    if let Some(recorder) = recorder.as_mut() {
        recorder.record(Event::new(
            "write_complete",
            None,
            json!({ "manifest": manifest.feature_records_path, "record_count": records.len() }),
        ))?;
        recorder.finish()?;
        println!("run_id={}", recorder.run_id());
    }

    println!("{}", serde_json::to_string_pretty(&manifest)?);
    Ok(())
}

fn cmd_graph(input: &Path, graph_type: &str, format: &str, out: Option<PathBuf>) -> Result<()> {
    if graph_type != "residue" {
        anyhow::bail!("only residue graph export is implemented in this bootstrap slice");
    }
    if format != "graphml" {
        anyhow::bail!("only graphml export is implemented in this bootstrap slice");
    }

    let structure = parse_file(input, ParseOptions::default())?;
    let graph = build_residue_graph(&structure, 8.0);
    let graphml = export_graphml(&graph);
    let out_path = normalize_output_path(out, &format!("{}.graphml", structure.id));
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    fs::write(&out_path, graphml).with_context(|| format!("failed to write {}", out_path.display()))?;
    println!("{}", out_path.display());
    Ok(())
}

fn cmd_provenance(sqlite: &Path) -> Result<()> {
    for line in inspect_sqlite(sqlite)? {
        println!("{line}");
    }
    Ok(())
}

fn collect_inputs(input: &Path) -> Result<Vec<PathBuf>> {
    if input.is_file() {
        return Ok(vec![input.to_path_buf()]);
    }

    let mut files = Vec::new();
    for entry in WalkDir::new(input) {
        let entry = entry?;
        if entry.file_type().is_file() && is_supported_file(entry.path()) {
            files.push(entry.path().to_path_buf());
        }
    }

    if files.is_empty() {
        anyhow::bail!("no supported structures found under {}", input.display());
    }

    files.sort();
    Ok(files)
}

fn is_supported_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| {
            let name = name.to_ascii_lowercase();
            name.ends_with(".pdb")
                || name.ends_with(".cif")
                || name.ends_with(".mmcif")
                || name.ends_with(".pdb.gz")
                || name.ends_with(".cif.gz")
                || name.ends_with(".mmcif.gz")
        })
        .unwrap_or(false)
}

#[allow(dead_code)]
fn _keep_structure_type(_: &Structure) {}
