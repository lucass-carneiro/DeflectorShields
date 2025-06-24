use crate::warp_drive_hamiltonian::WarpDriveHamiltonian;

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierreSharpRadii {
    radius: f64,
    sigma: f64,
}

#[derive(Debug, serde::Deserialize)]
pub enum WarpDriveAlcubierreSharpInteriorSpeed {
    Constant(f64),
    SmoothShutDown {
        initial_value: f64,
        shutdown_time: f64,
        shutdown_duration: f64,
    },
}

#[derive(Debug, serde::Deserialize)]
pub struct WarpDriveAlcubierreSharp {
    bubble_speed: f64,
    interior_speed: WarpDriveAlcubierreSharpInteriorSpeed,
    deflector_strength: f64,
    epsilon: f64,
    radii_x: WarpDriveAlcubierreSharpRadii,
    radii_y: WarpDriveAlcubierreSharpRadii,
}

impl WarpDriveAlcubierreSharp {
    pub fn r(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        f64::sqrt(
            self.epsilon
                + (-(t * self.bubble_speed) + x) * (-(t * self.bubble_speed) + x)
                + y * y
                + z * z,
        )
    }

    pub fn rho(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        y / f64::sqrt(self.epsilon + y * y + z * z)
    }

    pub fn f(&self, radii: &WarpDriveAlcubierreSharpRadii, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        if radii.radius < rr && (radii.radius + radii.sigma) > rr {
            1.0 / (1.0
                + f64::exp(
                    radii.sigma
                        * (1.0 / (radii.radius - rr) + 1.0 / (radii.radius - rr + radii.sigma)),
                ))
        } else if radii.radius < rr {
            0.0
        } else {
            1.0
        }
    }

    pub fn d_r_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        (self.bubble_speed * (t * self.bubble_speed - x))
            / f64::sqrt(
                self.epsilon
                    + (-(t * self.bubble_speed) + x) * (-(t * self.bubble_speed) + x)
                    + y * y
                    + z * z,
            )
    }

    pub fn d_r_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        (-(t * self.bubble_speed) + x)
            / f64::sqrt(
                self.epsilon
                    + (-(t * self.bubble_speed) + x) * (-(t * self.bubble_speed) + x)
                    + y * y
                    + z * z,
            )
    }

    pub fn d_r_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        y / f64::sqrt(
            self.epsilon
                + (-(t * self.bubble_speed) + x) * (-(t * self.bubble_speed) + x)
                + y * y
                + z * z,
        )
    }

    pub fn d_r_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (t, x, y, z) = (&q[0], &q[1], &q[2], &q[3]);
        z / f64::sqrt(
            self.epsilon
                + (-(t * self.bubble_speed) + x) * (-(t * self.bubble_speed) + x)
                + y * y
                + z * z,
        )
    }

    pub fn d_rho_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        (self.epsilon + z * z)
            / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z))
    }

    pub fn d_rho_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        let (y, z) = (&q[2], &q[3]);
        -((y * z) / ((self.epsilon + y * y + z * z) * f64::sqrt(self.epsilon + y * y + z * z)))
    }

    pub fn d_f_dr(&self, radii: &WarpDriveAlcubierreSharpRadii, q: &nalgebra::Vector4<f64>) -> f64 {
        let rr = self.r(q);

        if radii.radius < rr && (radii.radius + radii.sigma) > rr {
            -0.25
                * (radii.sigma
                    * (2.0 * ((radii.radius - rr) * (radii.radius - rr))
                        + 2.0 * (radii.radius - rr) * radii.sigma
                        + radii.sigma * radii.sigma)
                    * (1.0
                        / f64::cosh(
                            (radii.sigma
                                * (1.0 / (radii.radius - rr)
                                    + 1.0 / (radii.radius - rr + radii.sigma)))
                                / 2.,
                        )
                        * 1.0
                        / f64::cosh(
                            (radii.sigma
                                * (1.0 / (radii.radius - rr)
                                    + 1.0 / (radii.radius - rr + radii.sigma)))
                                / 2.,
                        )))
                / ((radii.radius - rr)
                    * (radii.radius - rr)
                    * ((radii.radius - rr + radii.sigma) * (radii.radius - rr + radii.sigma)))
        } else {
            0.0
        }
    }

    pub fn w(&self, t: f64) -> f64 {
        match self.interior_speed {
            WarpDriveAlcubierreSharpInteriorSpeed::Constant(speed) => speed,
            WarpDriveAlcubierreSharpInteriorSpeed::SmoothShutDown {
                initial_value,
                shutdown_time,
                shutdown_duration,
            } => {
                if t < shutdown_time {
                    initial_value
                } else if t > (shutdown_time + shutdown_duration) {
                    0.0
                } else {
                    (initial_value
                        * f64::powi(shutdown_duration + shutdown_time - t, 4)
                        * (f64::powi(shutdown_duration, 3)
                            + 10.0
                                * shutdown_duration
                                * ((shutdown_time - t) * (shutdown_time - t))
                            - 20.0 * f64::powi(shutdown_time - t, 3)
                            + 4.0 * (shutdown_duration * shutdown_duration) * (-shutdown_time + t)))
                        / f64::powi(shutdown_duration, 7)
                }
            }
        }
    }

    pub fn d_w_dt(&self, t: f64) -> f64 {
        match self.interior_speed {
            WarpDriveAlcubierreSharpInteriorSpeed::Constant(_) => 0.0,
            WarpDriveAlcubierreSharpInteriorSpeed::SmoothShutDown {
                initial_value,
                shutdown_time,
                shutdown_duration,
            } => {
                if shutdown_time > t || (shutdown_duration + shutdown_time) < t {
                    0.0
                } else {
                    (140.0
                        * initial_value
                        * f64::powi(shutdown_time - t, 3)
                        * f64::powi(shutdown_duration + shutdown_time - t, 3))
                        / f64::powi(shutdown_duration, 7)
                }
            }
        }
    }
}

impl WarpDriveHamiltonian for WarpDriveAlcubierreSharp {
    fn vx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.w(q[0]) * self.f(&self.radii_x, &q)
    }

    fn vy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.deflector_strength * self.rho(q) * self.f(&self.radii_y, &q)
    }

    fn vz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    // Vx Derivatives
    fn d_vx_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.w(q[0]) * self.d_f_dr(&self.radii_x, &q) * self.d_r_dt(&q)
            + self.d_w_dt(q[0]) * self.f(&self.radii_x, &q)
    }

    fn d_vx_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.w(q[0]) * self.d_f_dr(&self.radii_x, &q) * self.d_r_dx(&q)
    }

    fn d_vx_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.w(q[0]) * self.d_f_dr(&self.radii_x, &q) * self.d_r_dy(&q)
    }

    fn d_vx_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.w(q[0]) * self.d_f_dr(&self.radii_x, &q) * self.d_r_dz(&q)
    }

    // Vy Derivatives
    fn d_vy_dt(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.deflector_strength * self.rho(&q) * self.d_f_dr(&self.radii_y, &q) * self.d_r_dt(q)
    }

    fn d_vy_dx(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.deflector_strength * self.rho(&q) * self.d_f_dr(&self.radii_y, &q) * self.d_r_dx(q)
    }

    fn d_vy_dy(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.deflector_strength
            * (self.f(&self.radii_y, q) * self.d_rho_dy(q)
                + self.rho(q) * self.d_f_dr(&self.radii_y, q) * self.d_r_dy(q))
    }

    fn d_vy_dz(&self, q: &nalgebra::Vector4<f64>) -> f64 {
        self.deflector_strength
            * (self.f(&self.radii_y, q) * self.d_rho_dz(q)
                + self.rho(q) * self.d_f_dr(&self.radii_y, q) * self.d_r_dz(q))
    }

    // Vz derivatives
    fn d_vz_dt(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dx(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dy(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }

    fn d_vz_dz(&self, _: &nalgebra::Vector4<f64>) -> f64 {
        0.0
    }
}
