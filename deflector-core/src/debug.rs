use crate::types::ParticleState;
use crate::types::ParticleType;

pub fn is_normalized(
    state: &ParticleState<f64>,
    particle_type: &ParticleType,
    tolerance: f64,
) -> bool {
    let (vx, vy, vz, energy) = (state[3], state[4], state[5], state[6]);

    let energy_2 = energy * energy;
    let v_2 = vx * vx + vy * vy + vz * vz;
    let lhs = energy_2 * (1.0 - v_2);

    let eta = match particle_type {
        ParticleType::Massive => 1.0,
        ParticleType::Photon => 0.0,
    };

    f64::abs(lhs - eta) < tolerance
}
