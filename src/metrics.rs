use nalgebra::{Matrix4, Matrix4x2, Vector4};

pub trait Metric {
    fn g_ll(&self, x: Vector4<f64>) -> Matrix4<f64>;

    fn g_uu(&self, x: Vector4<f64>) -> Matrix4<f64> {
        self.g_ll(x)
            .try_inverse()
            .expect("Uninvertible spacetime metric")
    }

    fn hamiltonian(&self, p: Vector4<f64>, x: Vector4<f64>) -> f64 {
        let guu = self.g_uu(x);
        0.5 * (guu[(0, 0)] * p[0] * p[0]
            + 2.0 * guu[(0, 1)] * p[0] * p[1]
            + 2.0 * guu[(0, 2)] * p[0] * p[2]
            + 2.0 * guu[(0, 3)] * p[0] * p[3]
            + guu[(1, 1)] * p[1] * p[1]
            + 2.0 * guu[(1, 2)] * p[1] * p[2]
            + 2.0 * guu[(1, 3)] * p[1] * p[3]
            + guu[(2, 2)] * p[2] * p[2]
            + 2.0 * guu[(2, 3)] * p[2] * p[3]
            + guu[(3, 3)] * p[3] * p[3])
    }

    fn hamiltonian_derivs(&self, p: Vector4<f64>, x: Vector4<f64>) -> Matrix4x2<f64>;
}
