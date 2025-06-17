use std::fs::File;

use crate::types::ParticleState;
use polars::prelude::{Column, DataFrame};
use polars_io::SerWriter;
use polars_io::ipc::{IpcCompression, IpcWriter};

pub struct IpcMultiFile {
    df: DataFrame,
}

impl IpcMultiFile {
    pub fn new() -> Self {
        IpcMultiFile {
            df: DataFrame::default(),
        }
    }

    pub fn append(
        &mut self,
        particle_index: usize,
        iteration: usize,
        lambda: f64,
        state_vector: &ParticleState<f64>,
    ) {
        let col = Column::new(
            format!("particle {particle_index}").into(),
            [
                iteration as f64,
                lambda,
                state_vector[0],
                state_vector[1],
                state_vector[2],
                state_vector[3],
                state_vector[4],
                state_vector[5],
                state_vector[6],
                state_vector[7],
            ],
        );
        self.df.with_column(col).unwrap();
    }

    pub fn write(&mut self, file_name: &str, iteration: usize) {
        if !std::fs::exists(file_name).unwrap() {
            std::fs::create_dir(file_name).unwrap();
        }

        let mut file =
            File::create(format!("{file_name}/{file_name}_it{iteration:08}.ipc")).unwrap();

        IpcWriter::new(&mut file)
            .with_compression(Some(IpcCompression::LZ4))
            .with_parallel(true)
            .finish(&mut self.df)
            .unwrap();
    }
}
