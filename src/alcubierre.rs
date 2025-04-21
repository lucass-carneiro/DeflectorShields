use crate::metrics::Metric;
use nalgebra::{Matrix4, Matrix4x2, RowVector4, Vector4};

// https://arxiv.org/abs/gr-qc/9907019
#[derive(Debug)]
pub struct Alcubierre {
    v: f64,
    sigma: f64,
    r: f64,
}

impl Metric for Alcubierre {
    fn g_ll(&self, x: Vector4<f64>) -> Matrix4<f64> {
        let (t, x, y, z) = (x[0], x[1], x[2], x[3]);

        let r = f64::sqrt((x - self.v * t) + (x - self.v * t) + y * y + z * z);

        let f = (f64::tanh(self.sigma * (r + self.r)) - f64::tanh(self.sigma * (r - self.r)))
            / (2.0 * f64::tanh(self.sigma * self.r));

        let r1 = RowVector4::new(self.v * self.v * f * f - 1.0, -self.v * f, 0.0, 0.0);
        let r2 = RowVector4::new(-self.v * f, 1.0, 0.0, 0.0);
        let r3 = RowVector4::new(0.0, 0.0, 1.0, 0.0);
        let r4 = RowVector4::new(0.0, 0.0, 0.0, 1.0);

        Matrix4::from_rows(&[r1, r2, r3, r4])
    }

    fn hamiltonian_derivs(&self, p: Vector4<f64>, x: Vector4<f64>) -> Matrix4x2<f64> {
        let (t, x, y, z) = (x[0], x[1], x[2], x[3]);
        let (pt, px, py, pz) = (p[0], p[1], p[2], p[4]);

        let f = 0.0;
        let r = f64::sqrt((x - self.v * t) * (x - self.v * t) + y * y + z * z);

        let dh_dp = Vector4::new(-pt + px * self.v * f, p[1], p[2], p[3]);
        let dh_dq = Vector4::from_element(0.0);
        Matrix4x2::from_columns(&[dh_dp, dh_dq])
    }
}
