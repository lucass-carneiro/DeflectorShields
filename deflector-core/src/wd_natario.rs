use crate::errors::InitializationError;
use crate::types::{ParticleState, ParticleType};
use crate::warp_drive::WarpDrive;

#[derive(Debug, Clone)]
#[cfg_attr(feature = "jason_parfiles", derive(serde::Deserialize))]
pub struct WarpDriveNatario {
    pub radius: f64,
    pub sigma: f64,
    pub u: f64,
    pub x0: f64,
    pub t0: f64,
    pub epsilon: f64,
}

impl WarpDriveNatario {
    pub fn new(radius: f64, sigma: f64, u: f64) -> Self {
        WarpDriveNatario {
            radius,
            sigma,
            u,
            x0: 0.0,
            t0: 0.0,
            epsilon: 1.0e-12,
        }
    }

    pub fn resume(
        t: f64,
        radius: f64,
        sigma: f64,
        u: f64,
        x0: f64,
        t0: f64,
        epsilon: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Self {
        let mut wd = Self {
            radius,
            sigma,
            u,
            x0,
            t0,
            epsilon,
        };

        wd.adjust_ship_state(t, old_ship_state)
            .expect("Unable to resume warp drive state");

        wd
    }

    pub fn get_radius(&self) -> f64 {
        self.radius
    }

    pub fn get_sigma(&self) -> f64 {
        self.sigma
    }

    pub fn get_u(&self) -> f64 {
        self.u
    }

    pub fn update_radius(&mut self, new_radius: f64) {
        self.radius = new_radius;
    }

    pub fn update_sigma(&mut self, new_sigma: f64) {
        self.sigma = new_sigma;
    }

    pub fn update_u(&mut self, t: f64, new_u: f64) {
        self.x0 = self.get_bubble_position(t);
        self.t0 = t;
        self.u = new_u;
    }

    fn adjust_ship_state(
        &mut self,
        t: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Result<(), InitializationError> {
        self.x0 = old_ship_state[0];
        self.t0 = t;

        let new_ship_state = self.make_normalized_state(
            old_ship_state[0],
            0.0,
            0.0,
            self.u,
            0.0,
            0.0,
            &ParticleType::Massive,
        )?;

        *old_ship_state = new_ship_state;

        Ok(())
    }

    fn r(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }

    fn d_r_dx(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        (x - self.get_bubble_position(t))
            / f64::sqrt(
                self.epsilon
                    + f64::powi(y, 2)
                    + f64::powi(z, 2)
                    + f64::powi(x - self.get_bubble_position(t), 2),
            )
    }

    fn d_r_dy(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        y / f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }

    fn d_r_dz(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        z / f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }

    fn trans(&self, r: f64) -> f64 {
        if r < self.radius {
            1.0
        } else if r > (self.radius + self.sigma) {
            0.0
        } else {
            (f64::powi(-1. * r + self.radius + self.sigma, 4)
                * (20. * f64::powi(r - 1. * self.radius, 3)
                    + 10. * f64::powi(r - 1. * self.radius, 2) * self.sigma
                    + 4. * (r - 1. * self.radius) * f64::powi(self.sigma, 2)
                    + f64::powi(self.sigma, 3)))
                / f64::powi(self.sigma, 7)
        }
    }

    fn d_trans_dr(&self, r: f64) -> f64 {
        if r < self.radius || r > (self.radius + self.sigma) {
            0.0
        } else {
            (f64::powi(-1. * r + self.radius + self.sigma, 4)
                * (60. * f64::powi(r - 1. * self.radius, 2)
                    + 20. * (r - 1. * self.radius) * self.sigma
                    + 4. * f64::powi(self.sigma, 2)))
                / f64::powi(self.sigma, 7)
                - (4.
                    * f64::powi(-1. * r + self.radius + self.sigma, 3)
                    * (20. * f64::powi(r - 1. * self.radius, 3)
                        + 10. * f64::powi(r - 1. * self.radius, 2) * self.sigma
                        + 4. * (r - 1. * self.radius) * f64::powi(self.sigma, 2)
                        + f64::powi(self.sigma, 3)))
                    / f64::powi(self.sigma, 7)
        }
    }

    fn d2_trans_dr2(&self, r: f64) -> f64 {
        if r < self.radius || r > (self.radius + self.sigma) {
            0.0
        } else {
            (f64::powi(-1. * r + self.radius + self.sigma, 4)
                * (120. * (r - 1. * self.radius) + 20. * self.sigma))
                / f64::powi(self.sigma, 7)
                - (8.
                    * f64::powi(-1. * r + self.radius + self.sigma, 3)
                    * (60. * f64::powi(r - 1. * self.radius, 2)
                        + 20. * (r - 1. * self.radius) * self.sigma
                        + 4. * f64::powi(self.sigma, 2)))
                    / f64::powi(self.sigma, 7)
                + (12.
                    * f64::powi(-1. * r + self.radius + self.sigma, 2)
                    * (20. * f64::powi(r - 1. * self.radius, 3)
                        + 10. * f64::powi(r - 1. * self.radius, 2) * self.sigma
                        + 4. * (r - 1. * self.radius) * f64::powi(self.sigma, 2)
                        + f64::powi(self.sigma, 3)))
                    / f64::powi(self.sigma, 7)
        }
    }
}

impl WarpDrive for WarpDriveNatario {
    fn shut_down(&mut self, t: f64) {
        self.update_u(t, 0.0);
    }

    fn shut_up(
        &mut self,
        t: f64,
        old_ship_state: &mut crate::types::ParticleState<f64>,
        restart_parameters: &Self,
    ) -> Result<(), crate::errors::InitializationError> {
        self.radius = restart_parameters.get_radius();
        self.sigma = restart_parameters.get_sigma();
        self.u = restart_parameters.get_u();
        self.adjust_ship_state(t, old_ship_state)
    }

    fn get_bubble_position(&self, t: f64) -> f64 {
        self.x0 + self.u * (t - self.t0)
    }

    fn make_ship_state(
        &self,
        t: f64,
    ) -> Result<crate::types::ParticleState<f64>, crate::errors::InitializationError> {
        Ok(self.make_normalized_state(
            self.get_bubble_position(t),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            &ParticleType::Massive,
        )?)
    }

    fn update_ship_state(
        &self,
        t: f64,
        old_ship_state: &mut crate::types::ParticleState<f64>,
    ) -> Result<(), crate::errors::InitializationError> {
        let new_ship_state = self.make_ship_state(t)?;
        *old_ship_state = new_ship_state;
        Ok(())
    }

    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let lf = self.trans(lr);
        let dfdr = self.d_trans_dr(lr);
        (self.u * (2.0 * lf + (dfdr * (f64::powi(y, 2) + f64::powi(z, 2))) / lr)) / 2.
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        -0.5 * (dfdr * self.u * y * (x - self.get_bubble_position(t))) / lr
    }

    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        -0.5 * (dfdr * self.u * z * (x - self.get_bubble_position(t))) / lr
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdx = self.d_r_dx(t, x, y, z);
        (self.u
            * (2.0 * dfdr * drdx
                - (dfdr * drdx * (f64::powi(y, 2) + f64::powi(z, 2))) / f64::powi(lr, 2)
                + (d2fdr2 * drdx * (f64::powi(y, 2) + f64::powi(z, 2))) / lr))
            / 2.
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdy = self.d_r_dy(t, x, y, z);
        (self.u
            * (2.0 * dfdr * drdy + (2.0 * dfdr * y) / lr
                - (dfdr * drdy * (f64::powi(y, 2) + f64::powi(z, 2))) / f64::powi(lr, 2)
                + (d2fdr2 * drdy * (f64::powi(y, 2) + f64::powi(z, 2))) / lr))
            / 2.
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdz = self.d_r_dz(t, x, y, z);
        (self.u
            * (2.0 * dfdr * drdz + (2.0 * dfdr * z) / lr
                - (dfdr * drdz * (f64::powi(y, 2) + f64::powi(z, 2))) / f64::powi(lr, 2)
                + (d2fdr2 * drdz * (f64::powi(y, 2) + f64::powi(z, 2))) / lr))
            / 2.
    }

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdx = self.d_r_dx(t, x, y, z);
        -0.5 * (dfdr * self.u * y) / lr
            + (dfdr * drdx * self.u * y * (x - self.get_bubble_position(t)))
                / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdx * self.u * y * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdy = self.d_r_dy(t, x, y, z);
        -0.5 * (dfdr * self.u * (x - self.get_bubble_position(t))) / lr
            + (dfdr * drdy * self.u * y * (x - self.get_bubble_position(t)))
                / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdy * self.u * y * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdz = self.d_r_dz(t, x, y, z);
        (dfdr * drdz * self.u * y * (x - self.get_bubble_position(t))) / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdz * self.u * y * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdx = self.d_r_dx(t, x, y, z);
        -0.5 * (dfdr * self.u * z) / lr
            + (dfdr * drdx * self.u * z * (x - self.get_bubble_position(t)))
                / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdx * self.u * z * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdy = self.d_r_dy(t, x, y, z);
        (dfdr * drdy * self.u * z * (x - self.get_bubble_position(t))) / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdy * self.u * z * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let lr = self.r(t, x, y, z);
        let dfdr = self.d_trans_dr(lr);
        let d2fdr2 = self.d2_trans_dr2(lr);
        let drdz = self.d_r_dz(t, x, y, z);
        -0.5 * (dfdr * self.u * (x - self.get_bubble_position(t))) / lr
            + (dfdr * drdz * self.u * z * (x - self.get_bubble_position(t)))
                / (2. * f64::powi(lr, 2))
            - (d2fdr2 * drdz * self.u * z * (x - self.get_bubble_position(t))) / (2. * lr)
    }

    fn alp(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        1.0
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
