pub type ParticleState<T> = nalgebra::SVector<T, 8>;
pub type ParticleStates<T> = Vec<ParticleState<T>>;

#[derive(Clone, Copy, Debug, serde::Deserialize)]
pub enum ParticleType {
    Massive,
    Photon,
}
