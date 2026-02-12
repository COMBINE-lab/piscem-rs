use anyhow::Result;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedRadRecord {
    pub key: String,
}

pub trait RadDecoder {
    fn decode_records(&self, path: &str) -> Result<Vec<NormalizedRadRecord>>;
}

pub struct LibradiclAdapter;

impl RadDecoder for LibradiclAdapter {
    fn decode_records(&self, _path: &str) -> Result<Vec<NormalizedRadRecord>> {
        anyhow::bail!("libradicl adapter is a Phase 0 stub; wire develop-branch integration next")
    }
}
