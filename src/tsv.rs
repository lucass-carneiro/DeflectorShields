use std::fs::File;
use std::io::Write;

use crate::types::ParticleState;

#[derive(Debug)]
pub struct TSVFile {
    file: File,
}

impl TSVFile {
    pub fn new(file_name: &str) -> Self {
        let mut file = File::create(file_name).unwrap();
        writeln!(
            &mut file,
            "# 1: Iteration 2:Lambda 3:t 4:x 5:y 6:z 7:pt 8:px 9:py 10:pz"
        )
        .unwrap();
        TSVFile { file }
    }

    pub fn append(&mut self, iteration: u64, lambda: f64, state_vector: &ParticleState<f64>) {
        writeln!(
            self.file,
            "{}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}",
            iteration,
            lambda,
            state_vector[0],
            state_vector[1],
            state_vector[2],
            state_vector[3],
            state_vector[4],
            state_vector[5],
            state_vector[6],
            state_vector[7],
        )
        .unwrap();
    }
}
