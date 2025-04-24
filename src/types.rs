pub type ParticleState<T> = nalgebra::SVector<T, 8>;
//pub type PartivcleStates<T> = Vec<ParticleState<T>>;

#[derive(Debug, serde::Deserialize)]
pub enum ParticleType {
    Massive,
    Photon,
}
