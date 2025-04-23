use crate::metrics::Metric;
use nalgebra::{Matrix4, Matrix4x2, RowVector4, Vector4};

// https://arxiv.org/abs/gr-qc/9907019
#[derive(Debug)]
pub struct Alcubierre {
    v: f64,
    sigma: f64,
    radius: f64,
}

impl Alcubierre {
    pub fn new(v: f64, sigma: f64, radius: f64) -> Self {
        Alcubierre { v, sigma, radius }
    }

    pub fn r(&self, x: &Vector4<f64>) -> f64 {
        let (t, x, y, z) = (x[0], x[1], x[2], x[3]);
        f64::sqrt((x - self.v * t) * (x - self.v * t) + y * y + z * z)
    }

    pub fn f(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        (f64::tanh(self.sigma * (rr + self.radius)) - f64::tanh(self.sigma * (rr - self.radius)))
            / (2.0 * f64::tanh(self.sigma * self.radius))
    }

    pub fn dfdr(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        let a = 1.0 / f64::cosh((rr + self.radius) * self.sigma);
        let b = 1.0 / f64::cosh((rr - self.radius) * self.sigma);
        let c = f64::cosh(self.sigma * self.radius) / f64::sinh(self.sigma * self.radius);
        (self.sigma * c) * (a * a - b * b) / 2.0
    }

    pub fn drdt(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        let (t, x) = (x[0], x[1]);
        self.v * (t * self.v - x) / rr
    }

    pub fn drdx(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        let (t, x) = (x[0], x[1]);
        (x - self.v * t) / rr
    }

    pub fn drdy(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        x[2] / rr
    }

    pub fn drdz(&self, x: &Vector4<f64>) -> f64 {
        let rr = self.r(x);
        x[3] / rr
    }
}

impl Metric for Alcubierre {
    fn g_ll(&self, x: &Vector4<f64>) -> Matrix4<f64> {
        let ff = self.f(x);

        let r1 = RowVector4::new(self.v * self.v * ff * ff - 1.0, -self.v * ff, 0.0, 0.0);
        let r2 = RowVector4::new(-self.v * ff, 1.0, 0.0, 0.0);
        let r3 = RowVector4::new(0.0, 0.0, 1.0, 0.0);
        let r4 = RowVector4::new(0.0, 0.0, 0.0, 1.0);

        Matrix4::from_rows(&[r1, r2, r3, r4])
    }

    fn hamiltonian_derivs(&self, p: &Vector4<f64>, x: &Vector4<f64>) -> Matrix4x2<f64> {
        let (pt, px, py, pz) = (p[0], p[1], p[2], p[3]);

        let ff = self.f(x);
        let df_dr = self.dfdr(x);

        let dh_dpt = -pt - px * self.v * ff;
        let dh_dpx = px - self.v * ff * (pt + px * self.v * ff);
        let dh_dpy = py;
        let dh_dpz = pz;

        let dh_dt = -px * self.v * df_dr * self.drdt(x) * (pt + px * self.v * ff);
        let dh_dx = -px * self.v * df_dr * self.drdx(x) * (pt + px * self.v * ff);
        let dh_dy = -px * self.v * df_dr * self.drdy(x) * (pt + px * self.v * ff);
        let dh_dz = -px * self.v * df_dr * self.drdz(x) * (pt + px * self.v * ff);

        let dh_dp = Vector4::new(dh_dpt, dh_dpx, dh_dpy, dh_dpz);
        let dh_dq = Vector4::new(dh_dt, dh_dx, dh_dy, dh_dz);

        Matrix4x2::from_columns(&[dh_dp, dh_dq])
    }
}
