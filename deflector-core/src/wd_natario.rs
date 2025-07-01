use crate::transition::{cinf, d_cinf_dx, d2_cinf_dx2};
use crate::warp_drive::WarpDrive;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveNatario {
    pub radius: f64,
    pub sigma: f64,
    pub u: f64,
    pub epsilon: f64,
}

impl WarpDriveNatario {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let xi = x - self.u * t;
        f64::sqrt(xi * xi + y * y + z * z + self.epsilon)
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
}

impl WarpDrive for WarpDriveNatario {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, z) = (&q[0], &q[1], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);

        ((2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr)) * self.u * (t * self.u - x) * z) / (lr * lr)
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);

        (-dfdr - (2.0 * lf * (-1.0 + lr)) / (lr * lr)) * self.u * y * z
    }

    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let z = &q[3];
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);

        dfdr * self.u * (lr * lr - z * z)
            + (2.0 * lf * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z))) / (lr * lr)
    }

    // Vx Derivatives
    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, z) = (&q[0], &q[1], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        ((2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr)) * (self.u * self.u) * z) / (lr * lr)
            - (2.0
                * drdt
                * (2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / f64::powi(lr, 3)
            + ((2.0 * drdt * lf
                + 2.0 * dfdr * drdt * (-1.0 + lr)
                + 2.0 * dfdr * drdt * lr
                + d2fdr2 * drdt * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / (lr * lr)
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, z) = (&q[0], &q[1], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        -(((2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr)) * self.u * z) / (lr * lr))
            - (2.0
                * drdx
                * (2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / f64::powi(lr, 3)
            + ((2.0 * drdx * lf
                + 2.0 * dfdr * drdx * (-1.0 + lr)
                + 2.0 * dfdr * drdx * lr
                + d2fdr2 * drdx * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / (lr * lr)
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, z) = (&q[0], &q[1], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        (-2.0 * drdy * (2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr)) * self.u * (t * self.u - x) * z)
            / f64::powi(lr, 3)
            + ((2.0 * drdy * lf
                + 2.0 * dfdr * drdy * (-1.0 + lr)
                + 2.0 * dfdr * drdy * lr
                + d2fdr2 * drdy * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / (lr * lr)
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, z) = (&q[0], &q[1], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        ((2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr)) * self.u * (t * self.u - x)) / (lr * lr)
            - (2.0
                * drdz
                * (2.0 * lf * (-1.0 + lr) + dfdr * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / f64::powi(lr, 3)
            + ((2.0 * drdz * lf
                + 2.0 * dfdr * drdz * (-1.0 + lr)
                + 2.0 * dfdr * drdz * lr
                + d2fdr2 * drdz * (lr * lr))
                * self.u
                * (t * self.u - x)
                * z)
                / (lr * lr)
    }

    // Vy Derivatives
    fn d_vy_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        (-(d2fdr2 * drdt) + (4.0 * drdt * lf * (-1.0 + lr)) / f64::powi(lr, 3)
            - (2.0 * drdt * lf) / (lr * lr)
            - (2.0 * dfdr * drdt * (-1.0 + lr)) / (lr * lr))
            * self.u
            * y
            * z
    }

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        (-(d2fdr2 * drdx) + (4.0 * drdx * lf * (-1.0 + lr)) / f64::powi(lr, 3)
            - (2.0 * drdx * lf) / (lr * lr)
            - (2.0 * dfdr * drdx * (-1.0 + lr)) / (lr * lr))
            * self.u
            * y
            * z
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        (-dfdr - (2.0 * lf * (-1.0 + lr)) / (lr * lr)) * self.u * z
            + (-(d2fdr2 * drdy) + (4.0 * drdy * lf * (-1.0 + lr)) / f64::powi(lr, 3)
                - (2.0 * drdy * lf) / (lr * lr)
                - (2.0 * dfdr * drdy * (-1.0 + lr)) / (lr * lr))
                * self.u
                * y
                * z
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        (-dfdr - (2.0 * lf * (-1.0 + lr)) / (lr * lr)) * self.u * y
            + (-(d2fdr2 * drdz) + (4.0 * drdz * lf * (-1.0 + lr)) / f64::powi(lr, 3)
                - (2.0 * drdz * lf) / (lr * lr)
                - (2.0 * dfdr * drdz * (-1.0 + lr)) / (lr * lr))
                * self.u
                * y
                * z
    }

    // Vz derivatives
    fn d_vz_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let z = &q[3];
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdt = self.d_r_dt(&q);

        2.0 * dfdr * drdt * lr * self.u
            + d2fdr2 * drdt * self.u * (lr * lr - z * z)
            + (2.0 * lf * self.u * (3.0 * drdt * (lr * lr) - drdt * (z * z))) / (lr * lr)
            - (4.0 * drdt * lf * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z)))
                / f64::powi(lr, 3)
            + (2.0 * dfdr * drdt * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z))) / (lr * lr)
    }

    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let z = &q[3];
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        2.0 * dfdr * drdx * lr * self.u
            + d2fdr2 * drdx * self.u * (lr * lr - z * z)
            + (2.0 * lf * self.u * (3.0 * drdx * (lr * lr) - drdx * (z * z))) / (lr * lr)
            - (4.0 * drdx * lf * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z)))
                / f64::powi(lr, 3)
            + (2.0 * dfdr * drdx * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z))) / (lr * lr)
    }

    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let z = &q[3];
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        2.0 * dfdr * drdy * lr * self.u
            + d2fdr2 * drdy * self.u * (lr * lr - z * z)
            + (2.0 * lf * self.u * (3.0 * drdy * (lr * lr) - drdy * (z * z))) / (lr * lr)
            - (4.0 * drdy * lf * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z)))
                / f64::powi(lr, 3)
            + (2.0 * dfdr * drdy * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z))) / (lr * lr)
    }

    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let z = &q[3];
        let lr = self.r(&q);

        let lf = cinf(lr, 1.0, self.radius, self.sigma);
        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let d2fdr2 = d2_cinf_dx2(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        dfdr * self.u * (2.0 * drdz * lr - 2.0 * z)
            + d2fdr2 * drdz * self.u * (lr * lr - z * z)
            + (2.0
                * lf
                * self.u
                * (3.0 * drdz * (lr * lr) + 2.0 * z - 2.0 * lr * z - drdz * (z * z)))
                / (lr * lr)
            - (4.0 * drdz * lf * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z)))
                / f64::powi(lr, 3)
            + (2. * dfdr * drdz * self.u * (f64::powi(lr, 3) + z * z - lr * (z * z))) / (lr * lr)
    }

    fn alp(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        1.0
    }

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
