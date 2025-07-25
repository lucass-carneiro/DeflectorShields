use crate::dual::Dual;
use crate::errors::{InitializationError, SlippageError};
use crate::transition::*;
use crate::types::{ParticleState, ParticleType};
use crate::warp_drive::WarpDrive;

enum Transitions {
    Poly5,
    Poly7,
    Poly9,
    Cinf,
    Simple,
}

const TRANSITION: Transitions = Transitions::Poly7;

fn trans(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    match TRANSITION {
        Transitions::Poly5 => poly_trans_5(x.into(), y0.into(), x0.into(), dx.into()).f,
        Transitions::Poly7 => poly_trans_7(x.into(), y0.into(), x0.into(), dx.into()).f,
        Transitions::Poly9 => poly_trans_9(x.into(), y0.into(), x0.into(), dx.into()).f,
        Transitions::Cinf => cinf(x.into(), y0.into(), x0.into(), dx.into()).f,
        Transitions::Simple => simpl(x.into(), y0.into(), x0.into(), dx.into()).f,
    }
}

fn simpl(x: Dual, y0: Dual, x0: Dual, dx: Dual) -> Dual {
    let arg = (x - x0) / dx;
    y0 / (arg * arg + 1.0.into())
}

fn trans_dual(x: Dual, y0: f64, x0: f64, dx: f64) -> Dual {
    match TRANSITION {
        Transitions::Poly5 => poly_trans_5(x, y0.into(), x0.into(), dx.into()),
        Transitions::Poly7 => poly_trans_7(x, y0.into(), x0.into(), dx.into()),
        Transitions::Poly9 => poly_trans_9(x, y0.into(), x0.into(), dx.into()),
        Transitions::Cinf => cinf(x, y0.into(), x0.into(), dx.into()),
        Transitions::Simple => simpl(x, y0.into(), x0.into(), dx.into()),
    }
}

fn d_trans(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    match TRANSITION {
        Transitions::Poly5 => poly_trans_5(Dual::EPSILON + x.into(), y0.into(), x0.into(), dx.into()).df,
        Transitions::Poly7 => poly_trans_7(Dual::EPSILON + x.into(), y0.into(), x0.into(), dx.into()).df,
        Transitions::Poly9 => poly_trans_9(Dual::EPSILON + x.into(), y0.into(), x0.into(), dx.into()).df,
        Transitions::Cinf => cinf(Dual::EPSILON + x.into(), y0.into(), x0.into(), dx.into()).df,
        Transitions::Simple => simpl(Dual::EPSILON + x.into(), y0.into(), x0.into(), dx.into()).df,
    }
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "jason_parfiles", derive(serde::Deserialize))]
pub struct WarpDriveOurs {
    pub radius: f64,  // Bubble radius
    pub sigma: f64,   // Bubble transition region width
    pub u: f64,       // Bubble speed
    pub u0: f64,      // Dragging speed
    pub k0: f64,      // Deflection stength
    pub x0: f64,      // Initial bubble position
    pub t0: f64,      // Initial bubble time
    pub gamma: f64,   // Time dilation strength
    pub epsilon: f64, // Machine EPSILON
    pub deflector_sigma_pushout: f64,
    pub deflector_sigma_factor: f64,
    pub deflector_back: f64, // Should be 1.0 or 0.0
}

impl WarpDriveOurs {
    pub fn new(radius: f64, sigma: f64, u: f64, u0: f64, k0: f64, deflector_sigma_pushout: f64, deflector_sigma_factor: f64, deflector_back: f64) -> Self {
        WarpDriveOurs {
            radius,
            sigma,
            u,
            u0,
            k0,
            deflector_sigma_pushout: deflector_sigma_pushout,
            deflector_sigma_factor: deflector_sigma_factor,
            deflector_back: deflector_back,
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

    pub fn get_deflector_sigma_pushout(&self) -> f64 {
        self.deflector_sigma_pushout
    }

    pub fn get_deflector_sigma_factor(&self) -> f64 {
        self.deflector_sigma_factor
    }

    pub fn get_deflector_back(&self) -> f64 {
        self.deflector_back
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

        self.u0 += new_u - self.u;
        self.u = new_u;
    }

    pub fn update_u0(
        &mut self,
        new_u0: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Result<(), InitializationError> {
        self.u0 = new_u0;

        let vx = self.u - self.u0;

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

    pub fn update_deflector_sigma_pushout(&mut self, new_deflector_sigma_pushout: f64) {
        self.deflector_sigma_pushout = new_deflector_sigma_pushout;
    }

    pub fn update_deflector_sigma_factor(&mut self, new_deflector_sigma_factor: f64) {
        self.deflector_sigma_factor = new_deflector_sigma_factor;
    }

    pub fn update_deflector_back(&mut self, new_deflector_back: f64) {
        self.deflector_back = new_deflector_back;
    }

    fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }
    fn r_dual(&self, x: Dual, y: Dual, z: Dual) -> Dual {
        Dual::sqrt(
            self.epsilon
                + Dual::powi(y, 2)
                + Dual::powi(z, 2)
                + Dual::powi(x, 2),
        )
    }

    fn rho_y(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (q[2], q[3]);
        self.rho_y_dual(y.into(), z.into()).f
    }

    fn rho_y_dual(&self, y: Dual, z: Dual) -> Dual {
        y / Dual::sqrt(self.epsilon + Dual::powi(y, 2) + Dual::powi(z, 2))
    }

    fn rho_z(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (q[2], q[3]);
        z / f64::sqrt(self.epsilon + f64::powi(y, 2) + f64::powi(z, 2))
    }
    
    fn rho_z_dual(&self, y: Dual, z: Dual) -> Dual {
        z / Dual::sqrt(self.epsilon + Dual::powi(y, 2) + Dual::powi(z, 2))
    }

    fn d_r_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        let bubble_pos = self.get_bubble_position(t);
        (x - bubble_pos)
            / f64::sqrt(
            self.epsilon + f64::powi(y, 2) + f64::powi(z, 2) + f64::powi(x - bubble_pos, 2),
        )
    }

    fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        y / f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }

    fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (q[0], q[1], q[2], q[3]);
        z / f64::sqrt(
            self.epsilon
                + f64::powi(y, 2)
                + f64::powi(z, 2)
                + f64::powi(x - self.get_bubble_position(t), 2),
        )
    }

    fn d_rho_y_dy_old(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (q[2], q[3]);
        (self.epsilon + f64::powi(z, 2))
            / f64::powi(
            f64::sqrt(self.epsilon + f64::powi(y, 2) + f64::powi(z, 2)),
            3,
        )
    }
    fn d_rho_y_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.rho_y_dual(
            q[2] + Dual::EPSILON,
            q[3].into(),
        ).df
    }

    fn d_rho_y_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.rho_y_dual(
            q[2].into(),
            Dual::EPSILON + q[3].into(),
        ).df
    }

    fn d_rho_z_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (q[2], q[3]);
        -((y * z)
            / f64::powi(
            f64::sqrt(self.epsilon + f64::powi(y, 2) + f64::powi(z, 2)),
            3,
        ))
    }

    fn d_rho_z_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (q[2], q[3]);
        (self.epsilon + f64::powi(y, 2))
            / f64::powi(
            f64::sqrt(self.epsilon + f64::powi(y, 2) + f64::powi(z, 2)),
            3,
        )
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
        self.deflector_sigma_pushout = restart_parameters.get_deflector_sigma_pushout();
        self.deflector_sigma_factor = restart_parameters.get_deflector_sigma_factor();
        self.deflector_back = restart_parameters.get_deflector_back();

        // We force the position of the bubble to be where the ship is.
        // This is required because slippage can cause the ship to go back
        self.x0 = old_ship_state[0];
        self.t0 = t;

        // Compute the new ship state
        let vx = self.u - self.u0;

        if vx > 1.0 {
            return Err(InitializationError::from(SlippageError {
                u: self.u,
                u0: self.u0,
            }));
        }

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

        Ok(self.make_normalized_state(
            self.get_bubble_position(t),
            0.0,
            0.0,
            slip,
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
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let lr = self.r_dual(x.into(), y.into(), z.into());
        self.vx_dual(lr).f
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = y.into();
        let z_dual = z.into();

        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_y = self.rho_y_dual(y_dual, z_dual);

        self.vy_dual(x_dual, rho_y, lr).f
    }

    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = y.into();
        let z_dual = z.into();

        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_z = self.rho_z_dual(y_dual, z_dual);

        self.vz_dual(x_dual, rho_z, lr).f
    }
    // Vx Derivatives
    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let lr = self.r_dual(Dual::EPSILON + x.into(), y.into(), z.into());
        self.vx_dual(lr).df
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let lr = self.r_dual(x.into(), Dual::EPSILON + y.into(), z.into());
        self.vx_dual(lr).df
    }
    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let lr = self.r_dual(x.into(), y.into(), Dual::EPSILON + z.into());
        self.vx_dual(lr).df
    }

    // Vy Derivatives
    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = Dual::EPSILON + x.into();
        let y_dual = y.into();
        let z_dual = z.into();
        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_y = self.rho_y_dual(y_dual, z_dual);
        self.vy_dual(x_dual, rho_y, lr).df
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = Dual::EPSILON + y.into();
        let z_dual = z.into();
        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_y = self.rho_y_dual(y_dual, z_dual);
        self.vy_dual(x_dual, rho_y, lr).df
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = y.into();
        let z_dual = Dual::EPSILON + z.into();
        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_y = self.rho_y_dual(y_dual, z_dual);
        self.vy_dual(x_dual, rho_y, lr).df
    }

    // Vz derivatives
    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = Dual::EPSILON + x.into();
        let y_dual = y.into();
        let z_dual = z.into();

        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_z = self.rho_z_dual(y_dual, z_dual);

        self.vz_dual(x_dual, rho_z, lr).df
    }

    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = Dual::EPSILON + y.into();
        let z_dual = z.into();

        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_z = self.rho_z_dual(y_dual, z_dual);

        self.vz_dual(x_dual, rho_z, lr).df
    }

    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = y.into();
        let z_dual = Dual::EPSILON + z.into();

        let lr = self.r_dual(x_dual, y_dual, z_dual);
        let rho_z = self.rho_z_dual(y_dual, z_dual);

        self.vz_dual(x_dual, rho_z, lr).df
    }

    fn alp(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(q);
        let f = trans(lr, 1.0, self.radius, self.sigma);

        1.0 - self.gamma * f
    }

    fn d_alp_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(q);

        let dfdr = d_trans(lr, 1.0, self.radius, self.sigma);
        let drdx = self.d_r_dx(q);

        -self.gamma * dfdr * drdx
    }

    fn d_alp_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(q);

        let dfdr = d_trans(lr, 1.0, self.radius, self.sigma);
        let drdy = self.d_r_dy(q);

        -self.gamma * dfdr * drdy
    }

    fn d_alp_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let lr = self.r(q);

        let dfdr = d_trans(lr, 1.0, self.radius, self.sigma);
        let drdz = self.d_r_dz(q);

        -self.gamma * dfdr * drdz
    }

    fn vx_dual(&self, lr: Dual) -> Dual {
        let f = trans_dual(lr, 1.0, self.radius, self.sigma);
        self.u0 * f
    }

    fn v_base(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (x, y, z) = (q[1] - self.get_bubble_position(q[0]), q[2], q[3]);
        let x_dual = x.into();
        let y_dual = y.into();
        let z_dual = Dual::EPSILON + z.into();
        let lr = self.r_dual(x_dual, y_dual, z_dual);
        self.v_base_dual(x_dual, lr).f
    }

    fn v_base_dual(&self, x_dual: Dual, lr: Dual) -> Dual {
        let r0 = self.radius+self.deflector_sigma_pushout*self.sigma;
        let s0 = self.deflector_sigma_factor*self.sigma;
        let fa = trans_dual(lr - r0.into(), 1.0, 0.0, s0);
        let fb = trans_dual(-lr + r0.into(), 1.0, 0.0, s0);
        let fc = (1.0-self.deflector_back)*trans_dual(-x_dual+self.sigma.into(), 1.0, 0.0, self.sigma)+self.deflector_back.into();
        fa * fb* fc
    }

    fn vy_dual(&self, x_dual: Dual, rho_y: Dual, lr: Dual) -> Dual {
        self.k0 * rho_y * self.v_base_dual(x_dual, lr)
    }

    fn vz_dual(&self, x_dual: Dual, rho_z: Dual, lr: Dual) -> Dual {
        self.k0 * rho_z * self.v_base_dual(x_dual, lr)
    }
}
