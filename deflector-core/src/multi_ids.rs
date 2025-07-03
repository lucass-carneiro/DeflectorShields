use crate::errors::NormalizationError;
use crate::params::InitialData;
use crate::types::{ParticleStates, ParticleType};
use crate::warp_drive::WarpDrive;

use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

const RNG_SEED: [u8; 32] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
    27, 28, 29, 30, 31, 32,
];

fn make_single_particle<T: WarpDrive>(
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    particle_type: &ParticleType,
    warp_drive: &T,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut states: ParticleStates<f64> = Vec::new();

    let state = warp_drive.make_normalized_state(x, y, z, vx, vy, vz, particle_type)?;
    states.push(state);

    Ok(states)
}

fn make_random_y_stream<T: WarpDrive>(
    start: f64,
    length: f64,
    y_range: f64,
    vy_range: f64,
    num: u64,
    particle_type: &ParticleType,
    warp_drive: &T,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut rng = Xoshiro256PlusPlus::from_seed(RNG_SEED);

    let mut states: ParticleStates<f64> = Vec::new();

    let abs_length = f64::abs(length);
    let abs_height = f64::abs(y_range);
    let abs_vy = f64::abs(vy_range);

    let x0 = start;
    let xf = start + abs_length;
    let dx = (xf - x0) / (num as f64);

    let y0 = -abs_height;
    let yf = abs_height;

    let vy0 = -abs_vy;
    let vyf = abs_vy;

    // Particles
    for i in 0..num {
        let x = x0 + (i as f64) * dx;
        let y = rng.random_range(y0..=yf);
        let vy = rng.random_range(vy0..=vyf);

        let state = warp_drive.make_normalized_state(x, y, 0.0, 0.0, vy, 0.0, particle_type)?;
        states.push(state);
    }

    Ok(states)
}

pub fn make_initial_data<T: WarpDrive>(
    id: InitialData,
    particle_type: &ParticleType,
    warp_drive: &T,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut states = match id {
        InitialData::SingleParticle {
            x,
            y,
            z,
            vx,
            vy,
            vz,
        } => make_single_particle(x, y, z, vx, vy, vz, particle_type, warp_drive),
        InitialData::RandomYStream {
            start,
            length,
            y_range,
            vy_range,
            num,
        } => make_random_y_stream(
            start,
            length,
            y_range,
            vy_range,
            num,
            particle_type,
            warp_drive,
        ),
    }?;

    // Ship
    let ship_speed = warp_drive.get_bubble_speed() - warp_drive.get_dragging_speed();
    let ship = warp_drive
        .make_normalized_state(
            0.0,
            0.0,
            0.0,
            ship_speed / f64::sqrt(1.0 - ship_speed * ship_speed),
            0.0,
            0.0,
            &ParticleType::Massive,
        )
        .unwrap();
    states.push(ship);

    Ok(states)
}
