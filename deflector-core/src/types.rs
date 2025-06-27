pub type ParticleState<T> = nalgebra::SVector<T, 8>;
pub type ParticleStates<T> = Vec<ParticleState<T>>;

#[derive(Clone, Copy, Debug, serde::Deserialize)]
pub enum ParticleType {
    Massive,
    Photon,
}

pub trait ParticleStateComponents {
    type T;
    
    fn t(&self) -> Self::T;
    fn x(&self) -> Self::T;
    fn y(&self) -> Self::T;
    fn z(&self) -> Self::T;
    fn pt(&self) -> Self::T;
    fn px(&self) -> Self::T;
    fn py(&self) -> Self::T;
    fn pz(&self) -> Self::T;
}

impl <T: Copy> ParticleStateComponents for ParticleState<T> {
    type T = T;

    fn t(&self) -> Self::T {
        self[0]
    }

    fn x(&self) -> Self::T {
        self[1]
    }

    fn y(&self) -> Self::T {
        self[2]
    }

    fn z(&self) -> Self::T {
        self[3]
    }

    fn pt(&self) -> Self::T {
        self[4]
    }

    fn px(&self) -> Self::T {
        self[5]
    }

    fn py(&self) -> Self::T {
        self[6]
    }

    fn pz(&self) -> Self::T {
        self[7]
    }
}
