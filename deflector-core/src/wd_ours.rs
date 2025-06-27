use crate::transition::{cinf, d_cinf_dx};
use crate::warp_drive::WarpDrive;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveOurs {
    radius: f64,
    sigma: f64,
    u: f64,
    u0: f64,
    k0: f64,
    ts: f64,
    ds: f64,
    epsilon: f64,
}

impl WarpDriveOurs {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let xi = x - self.u * t;
        f64::sqrt(xi * xi + y * y + z * z + self.epsilon)
    }

    pub fn rho_y(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        y / f64::sqrt(y * y + z * z + self.epsilon)
    }

    pub fn rho_z(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        z / f64::sqrt(y * y + z * z + self.epsilon)
    }

    pub fn d_r_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        (self.u * (t * self.u - x))
            / f64::sqrt(self.epsilon + (-(t * self.u) + x) * (-(t * self.u) + x) + y * y + z * z)
    }

    pub fn d_r_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        (-(t * self.u) + x)
            / f64::sqrt(self.epsilon + (-(t * self.u) + x) * (-(t * self.u) + x) + y * y + z * z)
    }

    pub fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        y / f64::sqrt(self.epsilon + (-(t * self.u) + x) * (-(t * self.u) + x) + y * y + z * z)
    }

    pub fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        z / f64::sqrt(self.epsilon + (-(t * self.u) + x) * (-(t * self.u) + x) + y * y + z * z)
    }

    pub fn d_rho_y_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        (self.epsilon + z * z)
            / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z))
    }

    pub fn d_rho_y_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        -((y * z) / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z)))
    }

    pub fn d_rho_z_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        -((y * z) / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z)))
    }

    pub fn d_rho_z_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        (self.epsilon + y * y)
            / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z))
    }
}

impl WarpDrive for WarpDriveOurs {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let w = cinf(t, self.u0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        w * f
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let lrho = self.rho_y(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        w * lrho * f
    }

    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let lrho = self.rho_z(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        w * lrho * f
    }

    // Vx Derivatives
    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let w = cinf(t, self.u0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dwdt = d_cinf_dx(t, self.u0, self.ts, self.ds);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        dwdt * f + w * dfdr * drdt
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let w = cinf(t, self.u0, self.ts, self.ds);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        w * dfdr * drdx
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let w = cinf(t, self.u0, self.ts, self.ds);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        w * dfdr * drdy
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let w = cinf(t, self.u0, self.ts, self.ds);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        w * dfdr * drdz
    }

    // Vy Derivatives
    fn d_vy_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let rhoy = self.rho_y(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dwdt = d_cinf_dx(t, self.k0, self.ts, self.ds);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        dwdt * rhoy * f + w * rhoy * dfdr * drdt
    }

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let rhoy = self.rho_y(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        w * rhoy * dfdr * drdx
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let rhoy = self.rho_y(&q);
        let drhoydy = self.d_rho_y_dy(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        w * drhoydy * f + w * rhoy * dfdr * drdy
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let rhoy = self.rho_y(&q);
        let drhoydz = self.d_rho_y_dz(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        w * drhoydz * f + w * rhoy * dfdr * drdz
    }

    // Vz derivatives
    fn d_vz_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let rhoz = self.rho_z(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dwdt = d_cinf_dx(t, self.k0, self.ts, self.ds);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        dwdt * rhoz * f + w * rhoz * dfdr * drdt
    }

    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);
        let rhoz = self.rho_z(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        w * rhoz * dfdr * drdx
    }

    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let rhoz = self.rho_z(&q);
        let drhozdy = self.d_rho_z_dy(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        w * drhozdy * f + w * rhoz * dfdr * drdy
    }

    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let t = q[0];
        let lr = self.r(&q);

        let rhoz = self.rho_z(&q);
        let drhozdz = self.d_rho_z_dz(&q);

        let w = cinf(t, self.k0, self.ts, self.ds);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        w * drhozdz * f + w * rhoz * dfdr * drdz
    }

    // TODO: Create cool alpha
    fn alp(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        1.0
    }

    // TODO: Add cool alpha derivatives
    fn d_alp_dt(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_alp_dx(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_alp_dy(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_alp_dz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }
}
