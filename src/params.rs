use crate::errors::ParamError;

#[derive(Debug, serde::Deserialize)]
pub enum Normalization {
    Massive,
    Photon,
}

#[derive(Debug, serde::Deserialize)]
pub struct ParticleData {
    pub x1: f64,
    pub x2: f64,
    pub x3: f64,

    pub p1: f64,
    pub p2: f64,
    pub p3: f64,
}

#[derive(Debug, serde::Deserialize)]
pub struct TimeData {
    pub tf: f64,
    pub dt: f64,
}

#[derive(Debug, serde::Deserialize)]
pub struct SingleParams {
    pub particle: ParticleData,
    pub normalize_as: Normalization,
    pub time_integration: TimeData,
}

pub fn read_params(file_name: &str) -> Result<SingleParams, ParamError> {
    let par_file = match std::fs::read_to_string(file_name) {
        Ok(o) => o,
        Err(e) => {
            return Err(ParamError::FileIOError {
                file_name: String::from(file_name),
                io_error: e,
            });
        }
    };

    match serde_json::from_str::<SingleParams>(&par_file) {
        Ok(o) => Ok(o),
        Err(e) => Err(ParamError::FileParseError {
            file_name: String::from(file_name),
            json_error: e,
        }),
    }
}
