use crate::errors::ParamError;
use crate::types::ParticleType;

use crate::wd_ours::WarpDriveOurs;

#[derive(Debug, serde::Deserialize)]
pub struct AffineData {
    pub lambda_max: f64,
    pub dlambda: f64,
    pub out_every: Option<usize>,
}

#[derive(Debug, serde::Deserialize)]
pub enum MultiParticleID {
    StaticParticle {
        x: f64,
        y: f64,
        z: f64,
    },
    StaticWall {
        position: f64,
        extent: f64,
        num: u64,
    },
    StaticDebrisField {
        start: f64,
        width: f64,
        height: f64,
        num: u64,
    },
}

#[derive(Debug, serde::Deserialize)]
pub enum WarpDriveSolution {
    Ours(WarpDriveOurs),
}

#[derive(Debug, serde::Deserialize)]
pub struct Params {
    pub warp_drive_solution: WarpDriveSolution,
    pub normalize_as: ParticleType,
    pub affine_data: AffineData,
    pub multi_particle_id: MultiParticleID,
    pub ship_speed: f64,
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

    match serde_json::from_str::<Params>(&par_file) {
        Ok(o) => Ok(o),
        Err(e) => Err(ParamError::FileParseError {
            file_name: String::from(file_name),
            json_error: e,
        }),
    }
}
