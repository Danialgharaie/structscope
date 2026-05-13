use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuardConfig {
    pub enabled: bool,
}

pub fn guard_available() -> bool {
    false
}
