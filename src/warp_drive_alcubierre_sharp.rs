use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierreSharpRadii {
    radius: f64,
    sigma: f64,
}

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierreSharp {
    v: f64,
    k: f64,
    radii_x: WarpDriveAlcubierreSharpRadii,
    radii_y: WarpDriveAlcubierreSharpRadii,
}

impl WarpDriveAlcubierreSharp {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let p1 = x - self.v * t;
        f64::sqrt(p1 * p1 + y * y + z * z)
    }

    pub fn rho(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[1], &q[3]);
        y / f64::sqrt(y * y + z * z)
    }

    pub fn f(&self, radii: &WarpDriveAlcubierreSharpRadii, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);
        let radius = radii.radius;
        let sigma = radii.sigma;

        if radius < rr && rr < sigma {
            1.0 / (1.0
                + f64::exp(
                    ((radii.radius - radii.sigma) * (radii.radius - 2.0 * rr + radii.sigma))
                        / ((radii.radius - rr) * (rr - radii.sigma)),
                ))
        } else if radius < rr {
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

    pub fn d_rho_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let num = z * z;
        let den = f64::sqrt(y * y + z * z);
        num / (den * den * den)
    }

    pub fn d_rho_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let num = y * z;
        let den = f64::sqrt(y * y + z * z);
        -num / (den * den * den)
    }

    pub fn d_f_dr(&self, radii: &WarpDriveAlcubierreSharpRadii, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        if radii.radius < rr && rr < radii.sigma {
            (f64::exp(
                ((radii.radius - radii.sigma) * (radii.radius + radii.sigma))
                    / ((radii.radius - rr) * (rr - radii.sigma)),
            ) * (radii.radius - radii.sigma)
                * (radii.radius * radii.radius - 2.0 * radii.radius * rr + 2.0 * (rr * rr)
                    - 2.0 * rr * radii.sigma
                    + radii.sigma * radii.sigma))
                / ((f64::exp(
                    radii.radius / (rr - radii.sigma) + radii.sigma / (radii.radius - rr),
                ) + f64::exp(
                    radii.radius / (radii.radius - rr) + radii.sigma / (rr - radii.sigma),
                )) * (f64::exp(
                    radii.radius / (rr - radii.sigma) + radii.sigma / (radii.radius - rr),
                ) + f64::exp(
                    radii.radius / (radii.radius - rr) + radii.sigma / (rr - radii.sigma),
                )) * ((radii.radius - rr) * (radii.radius - rr))
                    * ((rr - radii.sigma) * (rr - radii.sigma)))
        } else {
            0.0
        }
    }
}

impl WarpDriveHamiltonian for WarpDriveAlcubierreSharp {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.v * self.f(&self.radii_x, &q)
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.k * (q[2] / self.rho(q)) * self.f(&self.radii_y, &q)
    }

    fn vz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    // Vx Derivatives
    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&self.radii_x, &q) * self.d_r_dt(&q)
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&self.radii_x, &q) * self.d_r_dx(&q)
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&self.radii_x, &q) * self.d_r_dy(&q)
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.d_f_dr(&self.radii_x, &q) * self.d_r_dz(&q)
    }

    // Vy Derivatives
    fn d_vy_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.k * (q[2] / self.rho(q)) * self.d_f_dr(&self.radii_y, &q) * self.d_r_dt(q)
    }

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.k * (q[2] / self.rho(q)) * self.d_f_dr(&self.radii_y, &q) * self.d_r_dx(q)
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.k
            * (self.f(&self.radii_y, q) * self.d_rho_dy(q)
                + self.rho(q) * self.d_f_dr(&self.radii_y, q) * self.d_r_dy(q))
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.k
            * (self.f(&self.radii_y, q) * self.d_rho_dz(q)
                + self.rho(q) * self.d_f_dr(&self.radii_y, q) * self.d_r_dz(q))
    }

    // Vz derivatives
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
