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

fn make_debris_field<T: WarpDrive>(
    x0: f64,
    width: f64,
    height: f64,
    vx_range: f64,
    vy_range: f64,
    num: u64,
    warp_drive: &T,
) -> Result<ParticleStates<f64>, NormalizationError> {
    let mut rng = Xoshiro256PlusPlus::from_seed(RNG_SEED);

    let mut states: ParticleStates<f64> = Vec::new();

    // We take abs of all parameters to make sure that even negative inputs
    // will produce correct results
    let abs_x0 = f64::abs(x0);
    let abs_width = f64::abs(width);
    let abs_height = f64::abs(height);
    let abs_vx = f64::abs(vx_range);
    let abs_vy = f64::abs(vy_range);

    let xf = abs_x0 + abs_width;

    let y0 = -abs_height / 2.0;
    let yf = abs_height / 2.0;

    let vx0 = -abs_vx;
    let vxf = abs_vx;

    let vy0 = -abs_vy;
    let vyf = abs_vy;

    // Particles
    for _ in 0..num {
        let x = rng.random_range(x0..=xf);
        let y = rng.random_range(y0..=yf);
        let vx = rng.random_range(vx0..=vxf);
        let vy = rng.random_range(vy0..=vyf);

        let state =
            warp_drive.make_normalized_state(x, y, 0.0, vx, vy, 0.0, &ParticleType::Massive)?;
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
        InitialData::DebrisField {
            x0,
            width,
            height,
            vx_range,
            vy_range,
            num,
        } => make_debris_field(x0, width, height, vx_range, vy_range, num, warp_drive),
    }?;

    // Ship
    states.push(warp_drive.make_ship_state(0.0).unwrap());

    Ok(states)
}
