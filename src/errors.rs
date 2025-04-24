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
#[error("State q = ({t}, {x}, {y}, {z}) p = (?, {px}, {py}, {pz}) is not normalizable")]
pub struct NormalizationError {
    pub particle_type: ParticleType,
    pub t: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub px: f64,
    pub py: f64,
    pub pz: f64,
}
