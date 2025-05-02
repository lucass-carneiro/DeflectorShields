use crate::types::ParticleState;
use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

pub fn rk4_step(
    h: f64,
    warp_drive_ham: &Box<dyn WarpDriveHamiltonian>,
    state: &mut ParticleState<f64>,
) {
    let y0 = state.clone();

    let k1 = warp_drive_ham.rhs(state);
    *state = y0 + h * k1 / 2.0;

    let k2 = warp_drive_ham.rhs(state);
    *state = y0 + h * k2 / 2.0;

    let k3 = warp_drive_ham.rhs(state);
    *state = y0 + h * k3;

    let k4 = warp_drive_ham.rhs(state);

    *state = y0 + h * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
}
