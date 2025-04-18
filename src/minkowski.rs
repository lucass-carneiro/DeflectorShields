use crate::metrics::Metric;
use nalgebra::{Matrix4, Matrix4x2, Vector4};

#[derive(Debug)]
pub struct Minkowski {}

impl Metric for Minkowski {
    fn g_ll(&self, _: Vector4<f64>) -> Matrix4<f64> {
        Matrix4::from_diagonal(&Vector4::new(-1.0, 1.0, 1.0, 1.0))
    }

    fn g_uu(&self, _: Vector4<f64>) -> Matrix4<f64> {
        Matrix4::from_diagonal(&Vector4::new(-1.0, 1.0, 1.0, 1.0))
    }

    fn hamiltonian(&self, p: Vector4<f64>, x: Vector4<f64>) -> f64 {
        let guu = self.g_uu(x);
        0.5 * (guu[(0, 0)] * p[0] * p[0]
            + 2.0 * guu[(0, 1)] * p[0] * p[1]
            + 2.0 * guu[(0, 2)] * p[0] * p[2]
            + 2.0 * guu[(0, 3)] * p[0] * p[3]
            + guu[(1, 1)] * p[1] * p[1]
            + 2.0 * guu[(1, 2)] * p[1] * p[2]
            + 2.0 * guu[(1, 3)] * p[1] * p[3]
            + guu[(2, 2)] * p[2] * p[2]
            + 2.0 * guu[(2, 3)] * p[2] * p[3]
            + guu[(3, 3)] * p[3] * p[3])
    }

    fn hamiltonian_derivs(&self, p: Vector4<f64>, x: Vector4<f64>) -> Matrix4x2<f64> {
        let dh_dp = Vector4::new(-p[0], p[1], p[2], p[3]);
        let dh_dq = Vector4::from_element(0.0);
        Matrix4x2::from_columns(&[dh_dp, dh_dq])
    }
}
