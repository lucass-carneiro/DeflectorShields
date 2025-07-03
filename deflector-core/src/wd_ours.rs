use crate::transition::{cinf, d_cinf_dx};
use crate::warp_drive::WarpDrive;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveOurs {
    pub radius: f64,
    pub sigma: f64,
    pub u: f64,
    pub u0: f64,
    pub k0: f64,
    pub ts: f64,
    pub ds: f64,
    pub gamma: f64,
    pub epsilon: f64,
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
    fn get_bubble_speed(&self) -> f64 {
        self.u
    }

    fn get_dragging_speed(&self) -> f64 {
        self.u0
    }

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

    fn alp(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);
        let f = cinf(lr, 1.0, self.radius, self.sigma);

        1.0 - self.gamma * f
    }

    fn d_alp_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        -self.gamma * dfdr * drdx
    }

    fn d_alp_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        -self.gamma * dfdr * drdy
    }

    fn d_alp_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        -self.gamma * dfdr * drdz
    }
}
