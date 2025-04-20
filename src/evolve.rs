use crate::errors::NormalizationError;
use crate::metrics::Metric;
use crate::params::{Normalization, SingleParams};

use nalgebra::Vector4;

#[derive(Debug)]
pub struct StateVector {
    pub p: Vector4<f64>,
    pub x: Vector4<f64>,
}

pub fn normalize<T: Metric>(
    params: SingleParams,
    metric: &T,
) -> Result<StateVector, NormalizationError> {
    let norm = match params.normalize_as {
        Normalization::Massive => -1.0,
        Normalization::Photon => 0.0,
    };

    let x = Vector4::new(
        0.0,
        params.particle.x1,
        params.particle.x2,
        params.particle.x3,
    );

    let guu = metric.g_uu(x);

    let mut p = Vector4::new(
        0.0,
        params.particle.p1,
        params.particle.p2,
        params.particle.p3,
    );

    let a = guu[(0, 0)];

    let b = 2.0 * (guu[(0, 1)] * p[1] + guu[(0, 2)] * p[2] + guu[(0, 3)] * p[3]);

    let c = guu[(1, 1)] * p[1] * p[1]
        + 2.0 * guu[(1, 2)] * p[1] * p[2]
        + 2.0 * guu[(1, 3)] * p[1] * p[3]
        + guu[(2, 2)] * p[2] * p[2]
        + 2.0 * guu[(2, 3)] * p[2] * p[3]
        + guu[(3, 3)] * p[3] * p[3]
        - norm;

    let delta = b * b - 4.0 * a * c;

    if delta < 0.0 {
        return Err(NormalizationError {
            normalization: params.normalize_as,
            p1: p[1],
            p2: p[2],
            p3: p[3],
        });
    }

    let p0_a = (-b - f64::sqrt(delta)) / (2.0 * a);
    let p0_b = (-b + f64::sqrt(delta)) / (2.0 * a);

    if p0_a > 0.0 {
        p[0] = p0_a;
    } else if p0_b > 0.0 {
        p[0] = p0_b;
    }

    Ok(StateVector { p, x })
}

pub fn rhs<T: Metric>(state: &StateVector, metric: &T) -> StateVector {
    let ham_derivs = metric.hamiltonian_derivs(state.p, state.x);

    let p_rhs = -Vector4::new(
        ham_derivs[(0, 1)],
        ham_derivs[(1, 1)],
        ham_derivs[(2, 1)],
        ham_derivs[(3, 1)],
    );

    let x_rhs = Vector4::new(
        ham_derivs[(0, 0)],
        ham_derivs[(1, 0)],
        ham_derivs[(2, 0)],
        ham_derivs[(3, 0)],
    );

    StateVector { p: p_rhs, x: x_rhs }
}

pub fn rk4_step<T: Metric>(state: &mut StateVector, metric: &T, h: f64) {
    let y0 = StateVector {
        p: state.p,
        x: state.x,
    };

    let k1 = rhs(state, metric);
    state.p = y0.p + h * k1.p / 2.0;
    state.x = y0.x + h * k1.x / 2.0;

    let k2 = rhs(state, metric);
    state.p = y0.p + h * k2.p / 2.0;
    state.x = y0.x + h * k2.x / 2.0;

    let k3 = rhs(state, metric);
    state.p = y0.p + h * k3.p;
    state.x = y0.x + h * k3.x;

    let k4 = rhs(state, metric);

    state.p = y0.p + h * (k1.p / 6.0 + k2.p / 3.0 + k3.p / 3.0 + k4.p / 6.0);
    state.x = y0.x + h * (k1.x / 6.0 + k2.x / 3.0 + k3.x / 3.0 + k4.x / 6.0);
}
