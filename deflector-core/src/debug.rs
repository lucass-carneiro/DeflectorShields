use crate::types::ParticleState;
use crate::types::ParticleType;

pub fn compute_norm(state: &ParticleState<f64>) -> f64 {
    let (vx, vy, vz, energy) = (state[3], state[4], state[5], state[6]);

    let energy_2 = energy * energy;
    let v_2 = vx * vx + vy * vy + vz * vz;
    energy_2 * (1.0 - v_2)
}

pub fn is_normalized(
    state: &ParticleState<f64>,
    particle_type: &ParticleType,
    tolerance: f64,
) -> bool {
    let lhs = compute_norm(state);

    let eta = match particle_type {
        ParticleType::Massive => 1.0,
        ParticleType::Photon => 0.0,
    };

    f64::abs(lhs - eta) < tolerance
}
