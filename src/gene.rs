pub struct Gene {
    pub start: u32,
    pub name: String,
}

impl Gene {
    pub fn get_name(&self) -> &str {
        &self.name
    }

    pub fn get_start(&self) -> u32 {
        self.start
    }
}
