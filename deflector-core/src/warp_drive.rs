use crate::errors::NormalizationError;
use crate::types::ParticleState;
use crate::types::ParticleType;

pub trait WarpDrive {
    fn get_bubble_speed(&self) -> f64;

    fn get_dragging_speed(&self) -> f64;

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
        let (vx_u, vy_u, vz_u) = (state[3], state[4], state[5]);
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

        let dx_udt = lvx + lalp * vx_u;

        let dy_udt = lvy + lalp * vy_u;

        let dz_udt = lvz + lalp * vz_u;

        let dvx_dt = (2.0 * dvxdx * vx_u
            + dvxdx * (lalp * lalp) * vx_u
            + dvxdy * lvx * vx_u
            + dvydx * lvx * vx_u
            + dvxdz * lvy * vx_u
            + dvzdx * lvy * vx_u
            + dvxdx * (lalp * lalp) * f64::powi(vx_u, 3)
            + dalpdx * (lalp * lalp) * (-1.0 + vx_u * vx_u)
            + dvxdy * vy_u
            + dvydx * vy_u
            + dvxdy * (lalp * lalp) * vy_u
            + 2.0 * dvydy * lvx * vy_u
            + dvydz * lvy * vy_u
            + dvzdy * lvy * vy_u
            + dalpdy * (lalp * lalp) * vx_u * vy_u
            + dvxdy * (lalp * lalp) * (vx_u * vx_u) * vy_u
            + dvydx * (lalp * lalp) * (vx_u * vx_u) * vy_u
            + dvydy * (lalp * lalp) * vx_u * (vy_u * vy_u)
            + (dvxdz
                + dvzdx
                + (dvydz + dvzdy) * lvx
                + 2.0 * dvzdz * lvy
                + dvxdz * (lalp * lalp) * (1.0 + vx_u * vx_u)
                + lalp * lalp * vx_u * (dalpdz + dvzdx * vx_u + (dvydz + dvzdy) * vy_u))
                * vz_u
            + dvzdz * (lalp * lalp) * vx_u * (vz_u * vz_u))
            / (lalp * lalp);

        let dvy_dt = -dalpdy - dvxdy * vx_u - dvydy * vy_u
            + dalpdx * vx_u * vy_u
            + dvxdx * (vx_u * vx_u) * vy_u
            + dalpdy * (vy_u * vy_u)
            + dvxdy * vx_u * (vy_u * vy_u)
            + dvydx * vx_u * (vy_u * vy_u)
            + dvydy * f64::powi(vy_u, 3)
            - dvzdy * vz_u
            + dalpdz * vy_u * vz_u
            + dvxdz * vx_u * vy_u * vz_u
            + dvzdx * vx_u * vy_u * vz_u
            + dvydz * (vy_u * vy_u) * vz_u
            + dvzdy * (vy_u * vy_u) * vz_u
            + dvzdz * vy_u * (vz_u * vz_u)
            + (lvx * lvx * ((dvxdy + dvydx) * vx_u + 2.0 * dvydy * vy_u + (dvydz + dvzdy) * vz_u))
                / (lalp * lalp)
            + (lvx
                * (2.0 * dvxdx * vx_u
                    + (dvxdz + dvzdx) * lvy * vx_u
                    + (dvxdy + dvydx) * vy_u
                    + (dvydz + dvzdy) * lvy * vy_u
                    + (dvxdz + dvzdx + 2.0 * dvzdz * lvy) * vz_u))
                / (lalp * lalp);

        let dvz_dt = -dalpdz
            + (lvy * lvy * ((dvxdz + dvzdx) * vx_u + (dvydz + dvzdy) * vy_u + 2.0 * dvzdz * vz_u))
                / (lalp * lalp)
            + dvxdz * vx_u * (-1.0 + vz_u * vz_u)
            + dvydz * vy_u * (-1.0 + vz_u * vz_u)
            + (lvy
                * (2.0 * dvxdx * vx_u
                    + (dvxdy + dvydx) * vy_u
                    + dvxdz * vz_u
                    + dvzdx * vz_u
                    + lvx
                        * ((dvxdy + dvydx) * vx_u + 2.0 * dvydy * vy_u + (dvydz + dvzdy) * vz_u)))
                / (lalp * lalp)
            + vz_u
                * (dalpdx * vx_u
                    + dvxdx * (vx_u * vx_u)
                    + vy_u * (dalpdy + (dvxdy + dvydx) * vx_u + dvydy * vy_u)
                    + (dalpdz + dvzdx * vx_u + dvzdy * vy_u) * vz_u
                    + dvzdz * (-1.0 + vz_u * vz_u));

        let de_dt = -(energy
            * (dalpdx * vx_u
                + dvxdx * (vx_u * vx_u)
                + vy_u * (dalpdy + (dvxdy + dvydx) * vx_u + dvydy * vy_u)
                + (dalpdz + (dvxdz + dvzdx) * vx_u + (dvydz + dvzdy) * vy_u) * vz_u
                + dvzdz * (vz_u * vz_u)));

        ParticleState::from_column_slice(&[dx_udt, dy_udt, dz_udt, dvx_dt, dvy_dt, dvz_dt, de_dt])
    }
}
