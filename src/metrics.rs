use nalgebra::{Matrix4, Matrix4x2, Vector4};

pub trait Metric {
    fn g_ll(&self, x: &Vector4<f64>) -> Matrix4<f64>;

    fn g_uu(&self, x: &Vector4<f64>) -> Matrix4<f64> {
        self.g_ll(x)
            .try_inverse()
            .expect("Uninvertible spacetime metric")
    }

    fn hamiltonian_derivs(&self, p: &Vector4<f64>, x: &Vector4<f64>) -> Matrix4x2<f64>;
}
