use std::fs::File;

use crate::types::ParticleState;
use polars::prelude::{DataFrame, df};
use polars_io::SerWriter;
use polars_io::ipc::{IpcCompression, IpcWriter};

pub struct IpcFile {
    df: DataFrame,
}

impl IpcFile {
    pub fn new() -> Self {
        let df = df!(
            "Iteration" => [0u64; 0],
            "Lambda" => [0.0; 0],
            "t" => [0.0; 0],
            "x" => [0.0; 0],
            "y" => [0.0; 0],
            "z" => [0.0; 0],
            "pt" => [0.0; 0],
            "px" => [0.0; 0],
            "py" => [0.0; 0],
            "pz" => [0.0; 0]
        )
        .unwrap();

        IpcFile { df }
    }

    pub fn append(&mut self, iteration: u64, lambda: f64, state_vector: &ParticleState<f64>) {
        let row_df = df!(
            "Iteration" => [iteration],
            "Lambda" => [lambda],
            "t" => [state_vector[0]],
            "x" => [state_vector[1]],
            "y" => [state_vector[2]],
            "z" => [state_vector[3]],
            "pt" => [state_vector[4]],
            "px" => [state_vector[5]],
            "py" => [state_vector[6]],
            "pz" => [state_vector[7]]
        )
        .unwrap();
        self.df.extend(&row_df).unwrap();
    }

    pub fn write(&mut self, file_name: &str) {
        let mut file = File::create(file_name).unwrap();
        IpcWriter::new(&mut file)
            .with_compression(Some(IpcCompression::LZ4))
            .with_parallel(true)
            .finish(&mut self.df)
            .unwrap();
    }
}
