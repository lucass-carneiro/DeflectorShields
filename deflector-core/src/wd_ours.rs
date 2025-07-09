use crate::errors::{InitializationError, SlippageError};
use crate::transition::{cinf, d_cinf_dx};
use crate::types::{ParticleState, ParticleType};
use crate::warp_drive::WarpDrive;

#[derive(Debug, Clone, serde::Deserialize)]
pub struct WarpDriveOurs {
    pub radius: f64,  // Bubble radius
    pub sigma: f64,   // Bubble transition region width
    pub u: f64,       // Bubble speed
    pub u0: f64,      // Dragging speed
    pub k0: f64,      // Deflection stength
    pub x0: f64,      // Initial bubble position
    pub t0: f64,      // Initial bubble time
    pub gamma: f64,   // Time dilation strength
    pub epsilon: f64, // Machine epsilon
}

impl WarpDriveOurs {
    pub fn new(radius: f64, sigma: f64, u: f64, u0: f64, k0: f64) -> Self {
        WarpDriveOurs {
            radius: radius,
            sigma: sigma,
            u: u,
            u0: u0,
            k0: k0,
            x0: 0.0,
            t0: 0.0,
            gamma: 0.0,
            epsilon: 1.0e-12,
        }
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

    pub fn get_u0(&self) -> f64 {
        self.u0
    }

    pub fn get_k0(&self) -> f64 {
        self.k0
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

        self.u0 = new_u - self.u + self.u0;
        self.u = new_u;
    }

    pub fn update_u0(
        &mut self,
        new_u0: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Result<(), InitializationError> {
        self.u0 = new_u0;

        let slip = self.u - self.u0;

        let vx = slip / f64::sqrt(1.0 - slip * slip);

        let new_ship_state = self.make_normalized_state(
            old_ship_state[0],
            0.0,
            0.0,
            vx,
            0.0,
            0.0,
            &ParticleType::Massive,
        )?;

        *old_ship_state = new_ship_state;

        Ok(())
    }

    pub fn update_k0(&mut self, new_k0: f64) {
        self.k0 = new_k0;
    }

    fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        f64::sqrt(
            self.epsilon
                + y * y
                + z * z
                + (x - self.get_bubble_position(t)) * (x - self.get_bubble_position(t)),
        )
    }

    fn rho_y(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        y / f64::sqrt(y * y + z * z + self.epsilon)
    }

    fn rho_z(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        z / f64::sqrt(y * y + z * z + self.epsilon)
    }

    fn d_r_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        (x - self.get_bubble_position(t))
            / f64::sqrt(
                self.epsilon
                    + y * y
                    + z * z
                    + (x - self.get_bubble_position(t)) * (x - self.get_bubble_position(t)),
            )
    }

    fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        y / f64::sqrt(
            self.epsilon
                + y * y
                + z * z
                + (x - self.get_bubble_position(t)) * (x - self.get_bubble_position(t)),
        )
    }

    fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        z / f64::sqrt(
            self.epsilon
                + y * y
                + z * z
                + (x - self.get_bubble_position(t)) * (x - self.get_bubble_position(t)),
        )
    }

    fn d_rho_y_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        (self.epsilon + z * z)
            / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z))
    }

    fn d_rho_y_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        -((y * z) / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z)))
    }

    fn d_rho_z_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        -((y * z) / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z)))
    }

    fn d_rho_z_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        (self.epsilon + y * y)
            / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z))
    }
}

impl WarpDrive for WarpDriveOurs {
    fn shut_down(&mut self, t: f64) {
        self.update_u(t, 0.0);

        // We don't call update_u0 because we want to preserve slippage
        self.u0 = 0.0;

        self.update_k0(0.0);
    }

    fn shut_up(
        &mut self,
        t: f64,
        old_ship_state: &mut ParticleState<f64>,
        restart_parameters: &Self,
    ) -> Result<(), InitializationError> {
        self.radius = restart_parameters.get_radius();
        self.sigma = restart_parameters.get_sigma();
        self.u = restart_parameters.get_u();
        self.u0 = restart_parameters.get_u0();
        self.k0 = restart_parameters.get_k0();

        // We force the position of the bubble to be where the ship is.
        // This is required because slippage can cause the ship to go back
        self.x0 = old_ship_state[0];
        self.t0 = t;

        // Compute the new ship state
        let slip = self.u - self.u0;

        if slip > 1.0 {
            return Err(InitializationError::from(SlippageError {
                u: self.u,
                u0: self.u0,
            }));
        }

        let vx = slip / f64::sqrt(1.0 - slip * slip);

        let new_ship_state = self.make_normalized_state(
            old_ship_state[0],
            0.0,
            0.0,
            vx,
            0.0,
            0.0,
            &ParticleType::Massive,
        )?;

        *old_ship_state = new_ship_state;

        Ok(())
    }

    fn get_bubble_position(&self, t: f64) -> f64 {
        self.x0 + self.u * (t - self.t0)
    }

    fn make_ship_state(&self, t: f64) -> Result<ParticleState<f64>, InitializationError> {
        let slip = self.u - self.u0;

        if slip > 1.0 {
            return Err(InitializationError::from(SlippageError {
                u: self.u,
                u0: self.u0,
            }));
        }

        let vx = slip / f64::sqrt(1.0 - slip * slip);

        Ok(self.make_normalized_state(
            self.get_bubble_position(t),
            0.0,
            0.0,
            vx,
            0.0,
            0.0,
            &ParticleType::Massive,
        )?)
    }

    fn update_ship_state(
        &self,
        t: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Result<(), InitializationError> {
        let new_ship_state = self.make_ship_state(t)?;
        *old_ship_state = new_ship_state;
        Ok(())
    }

    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        self.u0 * f
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);
        let lrho = self.rho_y(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        self.k0 * lrho * f
    }

    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);
        let lrho = self.rho_z(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        self.k0 * lrho * f
    }

    // Vx Derivatives
    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        self.u0 * dfdr * drdx
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        self.u0 * dfdr * drdy
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        self.u0 * dfdr * drdz
    }

    // Vy Derivatives
    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);
        let rhoy = self.rho_y(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        self.k0 * rhoy * dfdr * drdx
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let rhoy = self.rho_y(&q);
        let drhoydy = self.d_rho_y_dy(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        self.k0 * (drhoydy * f + rhoy * dfdr * drdy)
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let rhoy = self.rho_y(&q);
        let drhoydz = self.d_rho_y_dz(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        self.k0 * (drhoydz * f + rhoy * dfdr * drdz)
    }

    // Vz derivatives
    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);
        let rhoz = self.rho_z(&q);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(&q);

        self.k0 * rhoz * dfdr * drdx
    }

    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let rhoz = self.rho_z(&q);
        let drhozdy = self.d_rho_z_dy(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(&q);

        self.k0 * (drhozdy * f + rhoz * dfdr * drdy)
    }

    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(&q);

        let rhoz = self.rho_z(&q);
        let drhozdz = self.d_rho_z_dz(&q);

        let f = cinf(lr, 1.0, self.radius, self.sigma);

        let dfdr = d_cinf_dx(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(&q);

        self.k0 * (drhozdz * f + rhoz * dfdr * drdz)
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
