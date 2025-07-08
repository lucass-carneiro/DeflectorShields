use crate::errors::{InitializationError, SlippageError};
use crate::transition::{cinf, d_cinf_dx};
use crate::types::{ParticleState, ParticleType};
use crate::warp_drive::WarpDrive;

#[derive(Debug, Clone, serde::Deserialize)]
pub struct WarpDriveOurs {
    radius: f64,  // Bubble radius
    sigma: f64,   // Bubble transition region width
    u: f64,       // Bubble speed
    u0: f64,      // Dragging speed
    k0: f64,      // Deflection stength
    delta_u: f64, // TODO
    gamma: f64,   // Time dilation strength
    epsilon: f64, // Machine epsilon
}

impl WarpDriveOurs {
    pub fn new(radius: f64, sigma: f64, u: f64, u0: f64, k0: f64) -> Self {
        WarpDriveOurs {
            radius: radius,
            sigma: sigma,
            u: u,
            u0: u0,
            k0: k0,
            delta_u: 0.0,
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
        let udiff = new_u - self.u;
        self.delta_u = udiff * t + self.delta_u;
        self.u = new_u;
        self.u0 += udiff;
    }

    pub fn update_u0(
        &mut self,
        new_u0: f64,
        old_ship_state: &ParticleState<f64>,
    ) -> Result<ParticleState<f64>, InitializationError> {
        self.u0 = new_u0;

        let new_ship_speed = self.compute_ship_speed()?;

        let new_state = self.make_normalized_state(
            old_ship_state[0],
            old_ship_state[1],
            old_ship_state[2],
            new_ship_speed,
            old_ship_state[4],
            old_ship_state[5],
            &ParticleType::Massive,
        )?;

        Ok(new_state)
    }

    pub fn update_k0(&mut self, new_k0: f64) {
        self.k0 = new_k0;
    }
    
    pub fn shut_down_now(&mut self) {
        self.u0 = 0.;
        self.k0 = 0.;
        self.delta_u = 0.;
    }

    fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        f64::sqrt(
            self.epsilon
                + (self.delta_u - t * self.u + x) * (self.delta_u - t * self.u + x)
                + y * y
                + z * z,
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
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        (self.delta_u - t * self.u + x)
            / f64::sqrt(
                self.epsilon
                    + (self.delta_u - t * self.u + x) * (self.delta_u - t * self.u + x)
                    + y * y
                    + z * z,
            )
    }

    fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        y / f64::sqrt(
            self.epsilon
                + (self.delta_u - t * self.u + x) * (self.delta_u - t * self.u + x)
                + y * y
                + z * z,
        )
    }

    fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        z / f64::sqrt(
            self.epsilon
                + (self.delta_u - t * self.u + x) * (self.delta_u - t * self.u + x)
                + y * y
                + z * z,
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
    fn compute_ship_speed(&self) -> Result<f64, SlippageError> {
        let slip = self.u - self.u0;
        if slip > 1.0 {
            Err(SlippageError {
                u: self.u,
                u0: self.u0,
            })
        } else {
            Ok(slip / f64::sqrt(1.0 - slip * slip))
        }
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
