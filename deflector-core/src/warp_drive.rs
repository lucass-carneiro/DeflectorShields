use crate::errors::NormalizationError;
use crate::types::ParticleState;
use crate::types::ParticleType;

pub trait WarpDrive {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn vz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vy_dt(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_vz_dt(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vz_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vz_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_vz_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn alp(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn d_alp_dt(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_alp_dx(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_alp_dy(&self, q: &nalgebra::Vector4<f64>) -> f64;
    fn d_alp_dz(&self, q: &nalgebra::Vector4<f64>) -> f64;

    fn make_normalized_state(
        &self,
        t: f64,
        x: f64,
        y: f64,
        z: f64,
        px: f64,
        py: f64,
        pz: f64,
        particle_type: &ParticleType,
    ) -> Result<ParticleState<f64>, NormalizationError> {
        let q = nalgebra::Vector4::new(t, x, y, z);

        let lvx = self.vx(&q);
        let lvy = self.vy(&q);
        let lvz = self.vz(&q);

        let lalp = self.alp(&q);

        let delta = match particle_type {
            ParticleType::Massive => -1.0,
            ParticleType::Photon => 0.0,
        };

        let a = -(1.0 / (lalp * lalp));
        let b = (-2.0 * (lvx * px + lvy * py + lvz * pz)) / (lalp * lalp);
        let c = px * px + py * py + pz * pz
            - ((lvx * px + lvy * py + lvz * pz) * (lvx * px + lvy * py + lvz * pz)) / (lalp * lalp)
            - delta;

        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            Err(NormalizationError {
                particle_type: *particle_type,
                t,
                x,
                y,
                z,
                px,
                py,
                pz,
            })
        } else {
            let pt = (-b + f64::sqrt(discriminant)) / (2.0 * a);

            Ok(ParticleState::from_column_slice(&[
                t, x, y, z, pt, px, py, pz,
            ]))
        }
    }

    fn rhs(&self, state: &ParticleState<f64>) -> ParticleState<f64> {
        let q = nalgebra::Vector4::new(state[0], state[1], state[2], state[3]);
        let (pt, px, py, pz) = (&state[4], &state[5], state[6], state[7]);

        let lvx = self.vx(&q);
        let lvy = self.vy(&q);
        let lvz = self.vz(&q);

        let dvxdt = self.d_vx_dt(&q);
        let dvxdx = self.d_vx_dx(&q);
        let dvxdy = self.d_vx_dy(&q);
        let dvxdz = self.d_vx_dz(&q);

        let dvydt = self.d_vy_dt(&q);
        let dvydx = self.d_vy_dx(&q);
        let dvydy = self.d_vy_dy(&q);
        let dvydz = self.d_vy_dz(&q);

        let dvzdt = self.d_vz_dt(&q);
        let dvzdx = self.d_vz_dx(&q);
        let dvzdy = self.d_vz_dy(&q);
        let dvzdz = self.d_vz_dz(&q);

        let lalp = self.alp(&q);

        let dalpdt = self.d_alp_dt(&q);
        let dalpdx = self.d_alp_dx(&q);
        let dalpdy = self.d_alp_dy(&q);
        let dalpdz = self.d_alp_dz(&q);

        let dhdpt = -((pt + lvx * px + lvy * py + lvz * pz) / (lalp * lalp));

        let dhdpx = px - (lvx * (pt + lvx * px + lvy * py + lvz * pz)) / (lalp * lalp);

        let dhdpy = py - (lvy * (pt + lvx * px + lvy * py + lvz * pz)) / (lalp * lalp);

        let dhdpz = pz - (lvz * (pt + lvx * px + lvy * py + lvz * pz)) / (lalp * lalp);

        let dhdqt = ((pt + lvx * px + lvy * py + lvz * pz)
            * (-(lalp * (dvxdt * px + dvydt * py + dvzdt * pz))
                + dalpdt * (pt + lvx * px + lvy * py + lvz * pz)))
            / f64::powi(lalp, 3);

        let dhdqx = ((pt + lvx * px + lvy * py + lvz * pz)
            * (-(lalp * (dvxdx * px + dvydx * py + dvzdx * pz))
                + dalpdx * (pt + lvx * px + lvy * py + lvz * pz)))
            / f64::powi(lalp, 3);

        let dhdqy = ((pt + lvx * px + lvy * py + lvz * pz)
            * (-(lalp * (dvxdy * px + dvydy * py + dvzdy * pz))
                + dalpdy * (pt + lvx * px + lvy * py + lvz * pz)))
            / f64::powi(lalp, 3);

        let dhdqz = ((pt + lvx * px + lvy * py + lvz * pz)
            * (-(lalp * (dvxdz * px + dvydz * py + dvzdz * pz))
                + dalpdz * (pt + lvx * px + lvy * py + lvz * pz)))
            / f64::powi(lalp, 3);

        ParticleState::from_column_slice(&[
            dhdpt, dhdpx, dhdpy, dhdpz, -dhdqt, -dhdqx, -dhdqy, -dhdqz,
        ])
    }
}
