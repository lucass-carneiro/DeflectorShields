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

    fn hamiltonian_derivs(&self, p: Vector4<f64>, _: Vector4<f64>) -> Matrix4x2<f64> {
        let dh_dp = Vector4::new(-p[0], p[1], p[2], p[3]);
        let dh_dq = Vector4::from_element(0.0);
        Matrix4x2::from_columns(&[dh_dp, dh_dq])
    }
}
