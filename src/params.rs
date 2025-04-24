use crate::alcubierre::AlcubierreData;
use crate::errors::ParamError;
use crate::types::ParticleType;

#[derive(Debug, serde::Deserialize)]
pub struct AffineData {
    pub lambda_max: f64,
    pub dlambda: f64,
}

#[derive(Debug, serde::Deserialize)]
pub struct SingleParticleID {
    pub x0: f64,
    pub y0: f64,
    pub z0: f64,
    pub px0: f64,
    pub py0: f64,
    pub pz0: f64,
}

#[derive(Debug, serde::Deserialize)]
pub struct Params {
    pub alcubierre_data: AlcubierreData,
    pub normalize_as: ParticleType,
    pub affine_data: AffineData,
    pub single_particle_id: Option<SingleParticleID>,
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
