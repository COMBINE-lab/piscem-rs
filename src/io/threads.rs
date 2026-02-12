#[derive(Debug, Clone, Copy)]
pub struct ThreadConfig {
    pub threads: usize,
}

impl Default for ThreadConfig {
    fn default() -> Self {
        Self { threads: 1 }
    }
}
