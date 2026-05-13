use anyhow::Result;
use rusqlite::{params, Connection};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::fs;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use structscope_events::Event;
use uuid::Uuid;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvenanceConfig {
    pub sqlite_path: Option<PathBuf>,
    pub jsonl_path: Option<PathBuf>,
}

#[derive(Debug)]
pub struct ProvenanceRecorder {
    run_id: String,
    sqlite: Option<Connection>,
    jsonl_path: Option<PathBuf>,
}

impl ProvenanceRecorder {
    pub fn open(config: &ProvenanceConfig, command: &str) -> Result<Self> {
        let run_id = Uuid::new_v4().to_string();
        let sqlite = match &config.sqlite_path {
            Some(path) => {
                ensure_parent_dir(path)?;
                let conn = Connection::open(path)?;
                init_schema(&conn)?;
                conn.execute(
                    "INSERT INTO runs (run_id, command, started_at) VALUES (?1, ?2, datetime('now'))",
                    params![run_id, command],
                )?;
                Some(conn)
            }
            None => None,
        };

        Ok(Self {
            run_id,
            sqlite,
            jsonl_path: config.jsonl_path.clone(),
        })
    }

    pub fn run_id(&self) -> &str {
        &self.run_id
    }

    pub fn record(&mut self, event: Event) -> Result<()> {
        if let Some(conn) = &self.sqlite {
            conn.execute(
                "INSERT INTO events (run_id, event, timestamp, structure_id, details_json) VALUES (?1, ?2, ?3, ?4, ?5)",
                params![
                    &self.run_id,
                    &event.event,
                    event.timestamp.to_rfc3339(),
                    &event.structure_id,
                    event.details.to_string()
                ],
            )?;
        }

        if let Some(path) = &self.jsonl_path {
            append_jsonl(path, &event)?;
        }

        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.record(Event::new("run_complete", None, json!({ "run_id": &self.run_id })))
    }
}

pub fn inspect_sqlite(path: &Path) -> Result<Vec<String>> {
    let conn = Connection::open(path)?;
    let mut stmt = conn.prepare(
        "SELECT run_id, command, started_at FROM runs ORDER BY started_at DESC",
    )?;
    let rows = stmt.query_map([], |row| {
        Ok(format!(
            "{}\t{}\t{}",
            row.get::<_, String>(0)?,
            row.get::<_, String>(1)?,
            row.get::<_, String>(2)?
        ))
    })?;

    let mut out = Vec::new();
    for row in rows {
        out.push(row?);
    }
    Ok(out)
}

fn init_schema(conn: &Connection) -> Result<()> {
    conn.execute_batch(
        "
        CREATE TABLE IF NOT EXISTS runs (
            run_id TEXT PRIMARY KEY,
            command TEXT NOT NULL,
            started_at TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT NOT NULL,
            event TEXT NOT NULL,
            timestamp TEXT NOT NULL,
            structure_id TEXT,
            details_json TEXT NOT NULL
        );
        ",
    )?;
    Ok(())
}

fn append_jsonl(path: &Path, event: &Event) -> Result<()> {
    ensure_parent_dir(path)?;
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(file, "{}", serde_json::to_string(event)?)?;
    Ok(())
}

fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}
