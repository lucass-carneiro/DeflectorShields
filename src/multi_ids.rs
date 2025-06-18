use crate::errors::NormalizationError;
use crate::params::MultiParticleID;
use crate::types::{ParticleStates, ParticleType};
use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn make_static_wall(
    position: f64,
    extent: f64,
    num: u64,
    particle_type: &ParticleType,
    warp_drive_ham: &Box<dyn WarpDriveHamiltonian>,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut states: ParticleStates<f64> = Vec::new();

    let abs_extent = f64::abs(extent);
    let dy = (2.0 * abs_extent) / ((num - 1) as f64);

    // Particles
    for i in 0..num {
        let y = -abs_extent + (i as f64) * dy;
        let state = warp_drive_ham.make_normalized_state(
            0.0,
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

fn make_static_debris_field(
    start: f64,
    width: f64,
    height: f64,
    num: u64,
    particle_type: &ParticleType,
    warp_drive_ham: &Box<dyn WarpDriveHamiltonian>,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let seed: [u8; 32] = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
        26, 27, 28, 29, 30, 31, 32,
    ];

    let mut rng = Xoshiro256PlusPlus::from_seed(seed);

    let mut states: ParticleStates<f64> = Vec::new();

    let abs_width = f64::abs(width);
    let abs_height = f64::abs(height);

    let x0 = start;
    let xf = start + abs_width;

    let y0 = -abs_height / 2.0;
    let yf = abs_height / 2.0;

    // Particles
    for _ in 0..num {
        let x = rng.random_range(x0..=xf);
        let y = rng.random_range(y0..=yf);

        let state =
            warp_drive_ham.make_normalized_state(0.0, x, y, 0.0, 0.0, 0.0, 0.0, particle_type)?;
        states.push(state);
    }

    Ok(states)
}

pub fn make_multi_id(
    id: MultiParticleID,
    particle_type: &ParticleType,
    warp_drive_ham: &Box<dyn WarpDriveHamiltonian>,
) -> Result<ParticleStates<f64>, NormalizationError> {
    // Particles
    let mut states = match id {
        MultiParticleID::StaticWall {
            position,
            extent,
            num,
        } => make_static_wall(position, extent, num, particle_type, warp_drive_ham),

        MultiParticleID::StaticDebrisField {
            start,
            width,
            height,
            num,
        } => make_static_debris_field(start, width, height, num, particle_type, warp_drive_ham),
    }?;

    // Ship
    let ship = warp_drive_ham
        .make_normalized_state(
            0.0,
            0.0,
            0.0,
            0.0,
            warp_drive_ham.ship_speed(),
            0.0,
            0.0,
            particle_type,
        )
        .unwrap();
    states.push(ship);

    return Ok(states);
}
