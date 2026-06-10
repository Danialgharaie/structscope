use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use rayon::prelude::*;
use serde_json::json;
use std::fs;
use std::path::{Path, PathBuf};
use structscope_agent::guard_available;
use structscope_core::{kabsch, needleman_wunsch, parse_file, smith_waterman, three_to_one, ParseOptions, Structure};
use structscope_events::Event;
use structscope_features::{compute_features, per_residue::per_residue_features};
use structscope_graphs::{
    atom_id_to_residue_id, build_atom_graph, build_interface_graph, build_residue_graph,
    ChemicalInteraction, export_gml, export_graphml, export_html, export_json,
};
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
        #[arg(long, short = 'j')]
        jobs: Option<usize>,
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
    /// Optimal-superposition RMSD between two structures over matched atoms.
    Rmsd {
        reference: PathBuf,
        mobile: PathBuf,
        /// Atom selection for correspondence: ca, backbone, or all.
        #[arg(long, default_value = "ca")]
        atoms: String,
        /// Establish residue correspondence by sequence alignment (CA atoms);
        /// allows structures of different lengths.
        #[arg(long)]
        align: bool,
        /// Like --align but uses local (Smith-Waterman) alignment for partial
        /// or domain-level overlaps.
        #[arg(long)]
        local: bool,
    },
    /// Emit per-residue features (SASA, secondary structure, dihedrals) as JSONL.
    Residues {
        input: PathBuf,
        /// Optional output file; defaults to stdout.
        #[arg(long)]
        out: Option<PathBuf>,
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
            jobs,
        } => cmd_featurize(&input, &out, provenance, sqlite, jsonl, guard, jobs),
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
        Commands::Rmsd { reference, mobile, atoms, align, local } => cmd_rmsd(&reference, &mobile, &atoms, align, local),
        Commands::Residues { input, out } => cmd_residues(&input, out),
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
    jobs: Option<usize>,
) -> Result<()> {
    let inputs = collect_inputs(input)?;
    fs::create_dir_all(out)?;

    if let Some(num_threads) = jobs {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global();
    }

    if guard && !guard_available() {
        eprintln!("guard requested, but the optional eBPF agent is not implemented in this bootstrap slice");
    }

    let (tx, rx) = std::sync::mpsc::channel::<Event>();

    let logging_thread = if provenance {
        let mut recorder = ProvenanceRecorder::open(
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
        )?;

        let handle = std::thread::spawn(move || -> Result<ProvenanceRecorder> {
            while let Ok(event) = rx.recv() {
                recorder.record(event)?;
            }
            Ok(recorder)
        });
        Some(handle)
    } else {
        None
    };

    let tx_mutex = std::sync::Mutex::new(tx);

    let records: Vec<_> = inputs
        .par_iter()
        .filter_map(|path| {
            if provenance {
                let _ = tx_mutex.lock().unwrap().send(Event::new(
                    "parse_start",
                    None,
                    json!({ "path": path.display().to_string() }),
                ));
            }

            match parse_file(path, ParseOptions::default()) {
                Ok(structure) => {
                    let record = compute_features(&structure);
                    if provenance {
                        let _ = tx_mutex.lock().unwrap().send(Event::new(
                            "feature_complete",
                            Some(structure.id.clone()),
                            json!({ "path": path.display().to_string() }),
                        ));
                    }
                    Some(record)
                }
                Err(err) => {
                    if provenance {
                        let _ = tx_mutex.lock().unwrap().send(Event::new(
                            "structure_failed",
                            None,
                            json!({ "path": path.display().to_string(), "error": err.to_string() }),
                        ));
                    }
                    None
                }
            }
        })
        .collect();

    drop(tx_mutex);

    let mut recorder = if let Some(handle) = logging_thread {
        Some(handle.join().unwrap()?)
    } else {
        None
    };

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
    let ext = match format.to_lowercase().as_str() {
        "graphml" => "graphml",
        "gml" => "gml",
        "json" => "json",
        "html" => "html",
        other => anyhow::bail!("unknown format '{other}' (expected graphml, gml, json, or html)"),
    };

    let structure = parse_file(input, ParseOptions::default())?;

    let raw_interactions = structscope_features::interactions::interactions(&structure);
    let mut chemical_interactions = Vec::new();
    for ri in raw_interactions {
        if let (Some(res_a), Some(res_b)) = (atom_id_to_residue_id(&ri.atom_id_a), atom_id_to_residue_id(&ri.atom_id_b)) {
            chemical_interactions.push(ChemicalInteraction {
                kind: ri.kind.to_string(),
                res_id_a: res_a,
                res_id_b: res_b,
                distance: ri.distance,
            });
        }
    }

    let output = match (graph_type, ext) {
        ("residue", "graphml") => export_graphml(&build_residue_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("residue", "gml") => export_gml(&build_residue_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("residue", "json") => export_json(&build_residue_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("residue", "html") => {
            let pdb_data = structure_to_pdb(&structure);
            export_html(&build_residue_graph(&structure, 8.0, Some(&chemical_interactions)), &pdb_data, &structure.id)
        }

        ("interface", "graphml") => export_graphml(&build_interface_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("interface", "gml") => export_gml(&build_interface_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("interface", "json") => export_json(&build_interface_graph(&structure, 8.0, Some(&chemical_interactions))),
        ("interface", "html") => {
            let pdb_data = structure_to_pdb(&structure);
            export_html(&build_interface_graph(&structure, 8.0, Some(&chemical_interactions)), &pdb_data, &structure.id)
        }

        ("atom", "graphml") => export_graphml(&build_atom_graph(&structure, 5.0)),
        ("atom", "gml") => export_gml(&build_atom_graph(&structure, 5.0)),
        ("atom", "json") => export_json(&build_atom_graph(&structure, 5.0)),
        ("atom", "html") => {
            let pdb_data = structure_to_pdb(&structure);
            export_html(&build_atom_graph(&structure, 5.0), &pdb_data, &structure.id)
        }

        (other, _) => anyhow::bail!("unknown graph type '{other}' (expected residue, atom, or interface)"),
    };

    let out_path = normalize_output_path(out, &format!("{}.{}.{}", structure.id, graph_type, ext));
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    fs::write(&out_path, output).with_context(|| format!("failed to write {}", out_path.display()))?;
    println!("{}", out_path.display());
    Ok(())
}

fn structure_to_pdb(structure: &Structure) -> String {
    let mut out = String::new();
    let mut atom_serial = 1;
    for chain in &structure.chains {
        for residue in &chain.residues {
            for atom in &residue.atoms {
                let atom_name = if atom.name.len() < 4 {
                    format!(" {:<3}", atom.name)
                } else {
                    atom.name.clone()
                };
                let res_name = format!("{:<3}", residue.name);
                let chain_char = chain.id.chars().last().unwrap_or('A');
                let seq_num = residue.seq_number;
                let ins_code = residue.insertion_code.as_deref().unwrap_or(" ");
                let occupancy = atom.occupancy.unwrap_or(1.0);
                let temp_factor = atom.temp_factor.unwrap_or(0.0);
                let element = atom.element.as_deref().unwrap_or(" ");

                out.push_str(&format!(
                    "ATOM  {:5} {:4} {:3} {}{:4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:<2}\n",
                    atom_serial,
                    atom_name,
                    res_name,
                    chain_char,
                    seq_num,
                    ins_code,
                    atom.x,
                    atom.y,
                    atom.z,
                    occupancy,
                    temp_factor,
                    element
                ));
                atom_serial += 1;
            }
        }
    }
    out.push_str("END\n");
    out
}

fn cmd_rmsd(reference: &Path, mobile: &Path, atoms: &str, align: bool, local: bool) -> Result<()> {
    let reference_structure = parse_file(reference, ParseOptions::default())?;
    let mobile_structure = parse_file(mobile, ParseOptions::default())?;

    let (ref_coords, mob_coords) = if align || local {
        // Residue-level correspondence: align one-letter sequences, pair matched CA atoms.
        let residues = |s: &Structure| -> (Vec<u8>, Vec<[f64; 3]>) {
            let mut seq = Vec::new();
            let mut ca = Vec::new();
            for r in s.chains.iter().flat_map(|c| &c.residues) {
                if let Some(a) = r.atoms.iter().find(|a| a.name == "CA") {
                    seq.push(three_to_one(&r.name));
                    ca.push([a.x, a.y, a.z]);
                }
            }
            (seq, ca)
        };
        let (ref_seq, ref_ca) = residues(&reference_structure);
        let (mob_seq, mob_ca) = residues(&mobile_structure);
        let pairs = if local {
            smith_waterman(&ref_seq, &mob_seq)
        } else {
            needleman_wunsch(&ref_seq, &mob_seq)
        };
        let matched: Vec<(usize, usize)> = pairs.into_iter().filter(|&(i, j)| ref_seq[i] == mob_seq[j]).collect();
        if matched.is_empty() {
            anyhow::bail!("no matching residues found between the two structures");
        }
        (
            matched.iter().map(|&(i, _)| ref_ca[i]).collect::<Vec<_>>(),
            matched.iter().map(|&(_, j)| mob_ca[j]).collect::<Vec<_>>(),
        )
    } else {
        let select = |s: &Structure| -> Vec<[f64; 3]> {
            s.chains
                .iter()
                .flat_map(|c| &c.residues)
                .flat_map(|r| &r.atoms)
                .filter(|a| match atoms {
                    "ca" => a.name == "CA",
                    "backbone" => matches!(a.name.as_str(), "N" | "CA" | "C" | "O"),
                    _ => true,
                })
                .map(|a| [a.x, a.y, a.z])
                .collect()
        };
        let r = select(&reference_structure);
        let m = select(&mobile_structure);
        if r.len() != m.len() {
            anyhow::bail!(
                "atom count mismatch under selection '{atoms}': reference has {}, mobile has {} (use --align for sequence-based correspondence)",
                r.len(),
                m.len()
            );
        }
        (r, m)
    };

    let sp = kabsch(&mob_coords, &ref_coords).context("superposition failed (empty selection?)")?;
    let mode = if local {
        "local-ca"
    } else if align {
        "aligned-ca"
    } else {
        atoms
    };
    println!("rmsd={:.4}; atoms={}; selection={mode}", sp.rmsd, ref_coords.len());
    Ok(())
}

fn cmd_residues(input: &Path, out: Option<PathBuf>) -> Result<()> {
    let structure = parse_file(input, ParseOptions::default())?;
    let mut lines = String::new();
    for feature in per_residue_features(&structure) {
        lines.push_str(&serde_json::to_string(&feature)?);
        lines.push('\n');
    }
    match out {
        Some(path) => {
            fs::write(&path, lines).with_context(|| format!("failed to write {}", path.display()))?;
            println!("{}", path.display());
        }
        None => print!("{lines}"),
    }
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
                || name.ends_with(".bcif")
                || name.ends_with(".bcif.gz")
        })
        .unwrap_or(false)
}

#[allow(dead_code)]
fn _keep_structure_type(_: &Structure) {}
