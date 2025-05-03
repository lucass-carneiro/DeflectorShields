use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierreSharp {
    v: f64,
    sigma: f64,
    radius: f64,
}

impl WarpDriveAlcubierreSharp {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let p1 = x - self.v * t;
        f64::sqrt(p1 * p1 + y * y + z * z)
    }

    pub fn f(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        if self.radius < rr && rr < self.sigma {
            1.0 / (1.0
                + f64::exp(
                    ((self.radius - self.sigma) * (self.radius - 2.0 * rr + self.sigma))
                        / ((self.radius - rr) * (rr - self.sigma)),
                ))
        } else if self.radius < rr {
            0.0
        } else {
            1.0
        }
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

        if self.radius < rr && rr < self.sigma {
            let radius_minus_sigma = self.radius - self.sigma;
            let radius_minus_rr = self.radius - rr;
            let rr_minus_sigma = rr - self.sigma;
            let sech_term = 1.0
                / f64::cosh(
                    (radius_minus_sigma * (self.radius - 2.0 * rr + self.sigma))
                        / (2.0 * radius_minus_rr * rr_minus_sigma),
                );
            let radius_term = self.radius * self.radius - 2.0 * self.radius * rr + 2.0 * rr * rr
                - 2.0 * rr * self.sigma
                + self.sigma * self.sigma;
            (radius_minus_sigma * radius_term * sech_term * sech_term)
                / (4.0
                    * radius_minus_rr
                    * radius_minus_rr
                    * radius_minus_sigma
                    * radius_minus_sigma)
        } else {
            0.0
        }
    }
}

impl WarpDriveHamiltonian for WarpDriveAlcubierreSharp {
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
