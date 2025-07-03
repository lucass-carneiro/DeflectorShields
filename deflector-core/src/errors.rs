use crate::types::ParticleType;

#[derive(thiserror::Error, Debug)]
pub enum ParamError {
    #[error("Unable to read parameter file \"{file_name}\": {io_error}")]
    FileIOError {
        file_name: String,
        io_error: std::io::Error,
    },

    #[error("Unable to parse parameter file \"{file_name}\": {json_error}")]
    FileParseError {
        file_name: String,
        json_error: serde_json::Error,
    },
}

#[derive(thiserror::Error, Debug)]
#[error("State x = ({x}, {y}, {z}) v = (?, {vx}, {vy}, {vz}) is not normalizable")]
pub struct NormalizationError {
    pub particle_type: ParticleType,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}
