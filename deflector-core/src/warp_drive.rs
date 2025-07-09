use crate::errors::{InitializationError, NormalizationError};
use crate::types::ParticleState;
use crate::types::ParticleType;

pub trait WarpDrive {
    fn shut_down(&mut self, t: f64);

    fn shut_up(
        &mut self,
        t: f64,
        old_ship_state: &mut ParticleState<f64>,
        restart_parameters: Self,
    ) -> Result<(), InitializationError>;

    fn get_bubble_position(&self, t: f64) -> f64;

    fn make_ship_state(&self, t: f64) -> Result<ParticleState<f64>, InitializationError>;

    fn update_ship_state(
        &self,
        t: f64,
        old_ship_state: &mut ParticleState<f64>,
    ) -> Result<(), InitializationError>;

    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn alp(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_alp_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_alp_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_alp_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn make_normalized_state(
        &self,
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        particle_type: &ParticleType,
    ) -> Result<ParticleState<f64>, NormalizationError> {
        let delta = match particle_type {
            ParticleType::Massive => 1.0,
            ParticleType::Photon => 0.0,
        };

        let v2 = vx * vx + vy * vy + vz * vz;

        if v2 > 1.0 {
            return Err(NormalizationError {
                particle_type: *particle_type,
                x,
                y,
                z,
                vx,
                vy,
                vz,
            });
        }

        let energy = f64::sqrt(delta / (1.0 - v2));

        Ok(ParticleState::from_column_slice(&[
            x, y, z, vx, vy, vz, energy,
        ]))
    }

    fn rhs(&self, t: f64, state: &ParticleState<f64>) -> ParticleState<f64> {
        let q = nalgebra::Vector4::new(t, state[0], state[1], state[2]);
        let (velx, vely, velz) = (state[3], state[4], state[5]);
        let energy = state[6];

        let lvx = self.vx(&q);
        let lvy = self.vy(&q);
        let lvz = self.vz(&q);

        let dvxdx = self.d_vx_dx(&q);
        let dvxdy = self.d_vx_dy(&q);
        let dvxdz = self.d_vx_dz(&q);

        let dvydx = self.d_vy_dx(&q);
        let dvydy = self.d_vy_dy(&q);
        let dvydz = self.d_vy_dz(&q);

        let dvzdx = self.d_vz_dx(&q);
        let dvzdy = self.d_vz_dy(&q);
        let dvzdz = self.d_vz_dz(&q);

        let lalp = self.alp(&q);

        let dalpdx = self.d_alp_dx(&q);
        let dalpdy = self.d_alp_dy(&q);
        let dalpdz = self.d_alp_dz(&q);

        let dx_dt = lvx + lalp * velx;

        let dy_dt = lvy + lalp * vely;

        let dz_dt = lvz + lalp * velz;

        let dvelx_dt = dalpdx * (-1.0 + velx * velx) - dvydx * vely - dvzdx * velz
            + velx
                * (dvxdx * (-1.0 + velx * velx)
                    + vely * (dalpdy + (dvxdy + dvydx) * velx + dvydy * vely)
                    + (dalpdz + (dvxdz + dvzdx) * velx + (dvydz + dvzdy) * vely) * velz
                    + dvzdz * (velz * velz));

        let dvely_dt = dalpdy * (-1.0 + vely * vely) + dvxdy * velx * (-1.0 + vely * vely)
            - dvzdy * velz
            + vely
                * (dvydy * (-1.0 + vely * vely)
                    + dalpdz * velz
                    + (dvydz + dvzdy) * vely * velz
                    + dvzdz * (velz * velz)
                    + velx * (dalpdx + dvxdx * velx + dvydx * vely + (dvxdz + dvzdx) * velz));

        let dvelz_dt = -dalpdz - dvxdz * velx - dvydz * vely
            + (-dvzdz + velx * (dalpdx + dvxdx * velx)) * velz
            + vely * (dalpdy + (dvxdy + dvydx) * velx + dvydy * vely) * velz
            + (dalpdz + (dvxdz + dvzdx) * velx + (dvydz + dvzdy) * vely) * (velz * velz)
            + dvzdz * f64::powi(velz, 3);

        let denergy_dt = -(energy
            * (dalpdx * velx
                + dvxdx * (velx * velx)
                + vely * (dalpdy + (dvxdy + dvydx) * velx + dvydy * vely)
                + (dalpdz + (dvxdz + dvzdx) * velx + (dvydz + dvzdy) * vely) * velz
                + dvzdz * (velz * velz)));

        ParticleState::from_column_slice(&[
            dx_dt, dy_dt, dz_dt, dvelx_dt, dvely_dt, dvelz_dt, denergy_dt,
        ])
    }
}
