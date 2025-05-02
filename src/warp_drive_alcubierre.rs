use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierre {
    v: f64,
    sigma: f64,
    radius: f64,
}

impl WarpDriveAlcubierre {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let p1 = x - self.v * t;
        f64::sqrt(p1 * p1 + y * y + z * z)
    }

    pub fn f(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        (f64::tanh(self.sigma * (rr + self.radius)) - f64::tanh(self.sigma * (rr - self.radius)))
            / (2.0 * f64::tanh(self.sigma * self.radius))
    }

    pub fn d_r_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x) = (&q[0], &q[1]);
        let rr = self.r(q);

        self.v * (t * self.v - x) / rr
    }

    pub fn d_r_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x) = (&q[0], &q[1]);
        let rr = self.r(q);

        (x - self.v * t) / rr
    }

    pub fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);
        q[2] / rr
    }

    pub fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);
        q[3] / rr
    }

    pub fn d_f_dr(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        let coth_radius_sigma =
            f64::cosh(self.radius * self.sigma) / f64::sinh(self.radius * self.sigma);

        let sech_rr_p_radius = 1.0 / f64::cosh(self.sigma * (rr + self.radius));
        let sech_rr_m_radius = 1.0 / f64::cosh(self.sigma * (rr - self.radius));

        self.sigma
            * coth_radius_sigma
            * (sech_rr_p_radius * sech_rr_p_radius - sech_rr_m_radius * sech_rr_m_radius)
            / 2.0
    }
}

impl WarpDriveHamiltonian for WarpDriveAlcubierre {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.v * self.f(&q)
    }

    fn vy(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn vz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&q) * self.d_r_dt(&q)
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&q) * self.d_r_dx(&q)
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&q) * self.d_r_dy(&q)
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&q) * self.d_r_dz(&q)
    }

    fn d_vy_dt(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vy_dx(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vy_dy(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vy_dz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dt(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dx(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dy(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }
}
