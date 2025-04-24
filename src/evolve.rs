use crate::alcubierre::AlcubierreData;
use crate::types::ParticleState;

pub fn rk4_step(h: f64, alcubierre_data: &AlcubierreData, state: &mut ParticleState<f64>) {
    let y0 = state.clone();

    let k1 = alcubierre_data.rhs(state);
    *state = y0 + h * k1 / 2.0;

    let k2 = alcubierre_data.rhs(state);
    *state = y0 + h * k2 / 2.0;

    let k3 = alcubierre_data.rhs(state);
    *state = y0 + h * k3;

    let k4 = alcubierre_data.rhs(state);

    *state = y0 + h * (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0);
}
