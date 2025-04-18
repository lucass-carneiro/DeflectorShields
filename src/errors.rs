use crate::params::Normalization;

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
#[error("The 4 momentum ({p1}, {p2}, {p3}) is unnormalizable as a {normalization:?}")]
pub struct NormalizationError {
    pub normalization: Normalization,
    pub p1: f64,
    pub p2: f64,
    pub p3: f64,
}
