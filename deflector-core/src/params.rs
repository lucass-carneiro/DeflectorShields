use crate::errors::ParamError;
use crate::types::ParticleType;
use crate::wd_natario::WarpDriveNatario;
use crate::wd_ours::WarpDriveOurs;

#[derive(Debug, serde::Deserialize)]
pub struct TimeData {
    pub final_time: f64,
    pub time_step: f64,
    pub out_every: Option<usize>,
}

#[derive(Debug, serde::Deserialize)]
pub enum InitialData {
    SingleParticle {
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
    },
    DebrisField {
        x0: f64,
        width: f64,
        height: f64,
        vx_range: f64,
        vy_range: f64,
        num: u64,
    },
}

#[derive(Debug, serde::Deserialize)]
pub struct CommonParams {
    pub normalize_as: ParticleType,
    pub time_integration: TimeData,
    pub initial_data: InitialData,
}

#[derive(Debug, serde::Deserialize)]
pub enum WarpDriveType {
    Natario(WarpDriveNatario),
    Ours(WarpDriveOurs),
}

#[derive(Debug, serde::Deserialize)]
pub struct Params {
    pub warp_drive: WarpDriveType,
    pub normalize_as: ParticleType,
    pub time_integration: TimeData,
    pub initial_data: InitialData,
}

pub fn read_params(file_name: &str) -> Result<Params, ParamError> {
    let par_file = match std::fs::read_to_string(file_name) {
        Ok(o) => o,
        Err(e) => {
            return Err(ParamError::FileIOError {
                file_name: String::from(file_name),
                io_error: e,
            });
        }
    };

    match serde_json::from_str(&par_file) {
        Ok(o) => Ok(o),
        Err(e) => Err(ParamError::FileParseError {
            file_name: String::from(file_name),
            json_error: e,
        }),
    }
}
