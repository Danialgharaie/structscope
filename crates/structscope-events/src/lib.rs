use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Event {
    pub event: String,
    pub timestamp: DateTime<Utc>,
    pub structure_id: Option<String>,
    #[serde(default)]
    pub details: Value,
}

impl Event {
    pub fn new(event: impl Into<String>, structure_id: Option<String>, details: Value) -> Self {
        Self {
            event: event.into(),
            timestamp: Utc::now(),
            structure_id,
            details,
        }
    }
}
