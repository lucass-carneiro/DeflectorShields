use crate::errors::NormalizationError;
use crate::types::ParticleState;
use crate::types::ParticleType;

#[derive(Debug, serde::Deserialize)]
pub struct AlcubierreData {
    v: f64,
    sigma: f64,
    radius: f64,
}

impl AlcubierreData {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let p1 = x - self.v * t;
        f64::sqrt(p1 * p1 + y * y + z * z)
    }

    pub fn f(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        (f64::tanh(self.sigma * (rr + self.radius)) - f64::tanh(self.sigma * (rr - self.radius)))
            / (2.0 * f64::tanh(self.sigma * self.radius))
    }

    pub fn dr_dtxyz(&self, q: &nalgebra::Vector4<f64>) -> nalgebra::Vector4<f64> {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        let rr = self.r(q);

        nalgebra::Vector4::new(
            self.v * (t * self.v - x) / rr,
            (x - self.v * t) / rr,
            y / rr,
            z / rr,
        )
    }

    pub fn df_dr(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        let coth_radius_sigma =
            f64::cosh(self.radius * self.sigma) / f64::sinh(self.radius * self.sigma);

        let sech_rr_p_radius = 1.0 / f64::cosh(self.sigma * (rr + self.radius));
        let sech_rr_m_radius = 1.0 / f64::cosh(self.sigma * (rr - self.radius));

        self.sigma
            * coth_radius_sigma
            * (sech_rr_p_radius * sech_rr_p_radius - sech_rr_m_radius * sech_rr_m_radius)
            / 2.0
    }

    pub fn rhs(&self, state: &ParticleState<f64>) -> ParticleState<f64> {
        let q = nalgebra::Vector4::new(state[0], state[1], state[2], state[3]);
        let (pt, px, py, pz) = (&state[4], &state[5], state[6], state[7]);

        let fr = self.f(&q);
        let dfdr = self.df_dr(&q);
        let dtdtxyz = self.dr_dtxyz(&q);

        let dhdt = -(dfdr * dtdtxyz[0] * px * self.v * (pt + fr * px * self.v));
        let dhdx = -(dfdr * dtdtxyz[1] * px * self.v * (pt + fr * px * self.v));
        let dhdy = -(dfdr * dtdtxyz[2] * px * self.v * (pt + fr * px * self.v));
        let dhdz = -(dfdr * dtdtxyz[3] * px * self.v * (pt + fr * px * self.v));
        let dhdpt = -pt - fr * px * self.v;
        let dhdpx = px - fr * self.v * (pt + fr * px * self.v);
        let dhdpy = py;
        let dhdpz = pz;

        ParticleState::from_column_slice(&[dhdpt, dhdpx, dhdpy, dhdpz, -dhdt, -dhdx, -dhdy, -dhdz])
    }

    pub fn make_normalized_state(
        &self,
        x: f64,
        y: f64,
        z: f64,
        px: f64,
        py: f64,
        pz: f64,
        particle_type: ParticleType,
    ) -> Result<ParticleState<f64>, NormalizationError> {
        let fr = self.f(&nalgebra::Vector4::new(0.0, x, y, z));

        let delta = match particle_type {
            ParticleType::Massive => -1.0,
            ParticleType::Photon => 0.0,
        };

        let a = -0.5;
        let b = -fr * px * self.v;
        let c = (py * py + pz * pz + px * px * (1.0 - fr * fr * self.v * self.v) - delta) / 2.0;
        let delta = b * b - 4.0 * a * c;

        if delta < 0.0 {
            Err(NormalizationError {
                particle_type,
                t: 0.0,
                x,
                y,
                z,
                px,
                py,
                pz,
            })
        } else {
            let pt = (-b + f64::sqrt(delta)) / (2.0 * a);

            log::info!(
                "Normalized state q = ({}, {}, {}, {}) p = ({}, {}, {}, {})",
                0.0,
                x,
                y,
                z,
                pt,
                px,
                py,
                pz
            );

            Ok(ParticleState::from_column_slice(&[
                0.0, x, y, z, pt, px, py, pz,
            ]))
        }
    }
}
