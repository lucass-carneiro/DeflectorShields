use crate::debug::is_normalized;
use crate::types::{ParticleState, ParticleType};
use crate::warp_drive::WarpDrive;

pub fn rk4_step<T: WarpDrive>(
    time: f64,
    step_size: f64,
    warp_drive: &T,
    state: &mut ParticleState<f64>,
) {
    let y0 = *state;

    let k1 = warp_drive.rhs(time, state);
    *state = y0 + step_size * k1 / 2.0;

    let k2 = warp_drive.rhs(time + step_size / 2.0, state);
    *state = y0 + step_size * k2 / 2.0;

    let k3 = warp_drive.rhs(time + step_size / 2.0, state);
    *state = y0 + step_size * k3;

    let k4 = warp_drive.rhs(time + step_size, state);

    *state = y0 + step_size * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
}

pub fn rk4_step_adaptive<T: WarpDrive>(
    time: f64,
    step_size: f64,
    warp_drive: &T,
    state: &mut ParticleState<f64>,
    particle_type: &ParticleType
) {
    let y0 = *state;

    let k1 = warp_drive.rhs(time, state);
    *state = y0 + step_size * k1 / 2.0;

    let k2 = warp_drive.rhs(time + step_size / 2.0, state);
    *state = y0 + step_size * k2 / 2.0;

    let k3 = warp_drive.rhs(time + step_size / 2.0, state);
    *state = y0 + step_size * k3;

    let k4 = warp_drive.rhs(time + step_size, state);

    *state = y0 + step_size * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
    if step_size > 1e-8 && !is_normalized(state, particle_type,1e-6) {
        *state = y0;
        rk4_step(time, step_size / 2.0, warp_drive, state);
        rk4_step(time+step_size / 2.0, step_size / 2.0, warp_drive, state);
    }
}
