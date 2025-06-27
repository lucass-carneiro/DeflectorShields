use crate::types::ParticleState;
use crate::warp_drive::WarpDrive;

pub fn rk4_step_dyn(step_size: f64, warp_drive: &Box<dyn WarpDrive>, state: &mut ParticleState<f64>) {
    let y0 = *state;

    let k1 = warp_drive.rhs(state);
    *state = y0 + step_size * k1 / 2.0;

    let k2 = warp_drive.rhs(state);
    *state = y0 + step_size * k2 / 2.0;

    let k3 = warp_drive.rhs(state);
    *state = y0 + step_size * k3;

    let k4 = warp_drive.rhs(state);

    *state = y0 + step_size * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
}

pub fn rk4_step<T: WarpDrive>(step_size: f64, warp_drive: &T, state: &mut ParticleState<f64>) {
    let y0 = *state;

    let k1 = warp_drive.rhs(state);
    *state = y0 + step_size * k1 / 2.0;

    let k2 = warp_drive.rhs(state);
    *state = y0 + step_size * k2 / 2.0;

    let k3 = warp_drive.rhs(state);
    *state = y0 + step_size * k3;

    let k4 = warp_drive.rhs(state);

    *state = y0 + step_size * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
}
