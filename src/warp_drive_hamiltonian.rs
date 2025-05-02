use crate::errors::NormalizationError;
use crate::types::ParticleState;
use crate::types::ParticleType;

pub trait WarpDriveHamiltonian {
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

    // Position is contra-variant. Momentum is co-variant
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

        let delta = match particle_type {
            ParticleType::Massive => -1.0,
            ParticleType::Photon => 0.0,
        };

        let a = -0.5;
        let b = -(lvx * px + lvy * py + lvz * pz);
        let c = (-((-1.0 + lvx * lvx) * px * px)
            - (-1.0 + lvy * lvy) * py * py
            - 2.0 * lvy * lvz * py * pz
            + pz * pz
            - lvz * lvz * pz * pz * -2.0 * lvx * px * (lvy * py + lvz * pz)
            - delta)
            / 2.0;
        let discriminant = b * b - 4.0 * a * c;

        if discriminant < 0.0 {
            Err(NormalizationError {
                particle_type: *particle_type,
                t: 0.0,
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
                0.0, x, y, z, pt, px, py, pz,
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

        let dhdpt = -pt - lvx * px - lvy * py - lvz * pz;
        let dhdpx = px - lvx * (pt + lvx * px + lvy * py + lvz * pz);
        let dhdpy = py - lvy * (pt + lvx * px + lvy * py + lvz * pz);
        let dhdpz = pz - lvz * (pt + lvx * px + lvy * py + lvz * pz);

        let dhdqt =
            -((dvxdt * px + dvydt * py + dvzdt * pz) * (pt + lvx * px + lvy * py + lvz * pz));

        let dhdqx =
            -((dvxdx * px + dvydx * py + dvzdx * pz) * (pt + lvx * px + lvy * py + lvz * pz));

        let dhdqy =
            -((dvxdy * px + dvydy * py + dvzdy * pz) * (pt + lvx * px + lvy * py + lvz * pz));

        let dhdqz =
            -((dvxdz * px + dvydz * py + dvzdz * pz) * (pt + lvx * px + lvy * py + lvz * pz));

        ParticleState::from_column_slice(&[
            dhdpt, dhdpx, dhdpy, dhdpz, -dhdqt, -dhdqx, -dhdqy, -dhdqz,
        ])
    }
}
