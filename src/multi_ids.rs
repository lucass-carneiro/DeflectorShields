use crate::alcubierre::AlcubierreData;
use crate::errors::NormalizationError;
use crate::params::MultiParticleID;
use crate::types::{ParticleStates, ParticleType};

fn make_static_wall(
    position: f64,
    extent: f64,
    num: u64,
    particle_type: &ParticleType,
    alcubierre_data: &AlcubierreData,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut states: ParticleStates<f64> = Vec::new();

    let abs_extent = f64::abs(extent);
    let dy = (2.0 * abs_extent) / ((num - 1) as f64);

    for i in 0..num {
        let y = -abs_extent + (i as f64) * dy;
        let state = alcubierre_data.make_normalized_state(
            position,
            y,
            0.0,
            0.0,
            0.0,
            0.0,
            particle_type,
        )?;
        states.push(state);
    }

    Ok(states)
}

pub fn make_multi_id(
    id: MultiParticleID,
    particle_type: &ParticleType,
    alcubierre_data: &AlcubierreData,
) -> Result<ParticleStates<f64>, NormalizationError> {
    match id {
        MultiParticleID::StaticWall {
            position,
            extent,
            num,
        } => make_static_wall(position, extent, num, particle_type, alcubierre_data),
    }
}
