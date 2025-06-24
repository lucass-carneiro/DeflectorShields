pub type ParticleState<T> = nalgebra::SVector<T, 8>;
pub type ParticleStates<T> = Vec<ParticleState<T>>;

#[derive(Clone, Copy, Debug, serde::Deserialize)]
pub enum ParticleType {
    Massive,
    Photon,
}

#[derive(Clone, Copy, Debug)]
pub struct ParticleStateComponents<T> {
    pub t: T,
    pub x: T,
    pub y: T,
    pub z: T,
    pub pt: T,
    pub px: T,
    pub py: T,
    pub pz: T,
}

impl<T: Copy> ParticleStateComponents<T> {
    pub fn from_state(state: &ParticleState<T>) -> ParticleStateComponents<T> {
        ParticleStateComponents {
            t: state[0],
            x: state[0],
            y: state[0],
            z: state[0],
            pt: state[0],
            px: state[0],
            py: state[0],
            pz: state[0],
        }
    }
}
