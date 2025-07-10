pub type ParticleState<T> = nalgebra::SVector<T, 7>;
pub type ParticleStates<T> = Vec<ParticleState<T>>;

#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "jason_parfiles", derive(serde::Deserialize))]
pub enum ParticleType {
    Massive,
    Photon,
}

pub trait ParticleStateComponents {
    type T;

    fn x(&self) -> Self::T;
    fn y(&self) -> Self::T;
    fn z(&self) -> Self::T;
    fn vx(&self) -> Self::T;
    fn vy(&self) -> Self::T;
    fn vz(&self) -> Self::T;
    fn energy(&self) -> Self::T;
}

impl<T: Copy> ParticleStateComponents for ParticleState<T> {
    type T = T;

    fn x(&self) -> Self::T {
        self[0]
    }

    fn y(&self) -> Self::T {
        self[1]
    }

    fn z(&self) -> Self::T {
        self[2]
    }

    fn vx(&self) -> Self::T {
        self[3]
    }

    fn vy(&self) -> Self::T {
        self[4]
    }

    fn vz(&self) -> Self::T {
        self[5]
    }

    fn energy(&self) -> Self::T {
        self[6]
    }
}
