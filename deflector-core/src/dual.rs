use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Copy, Clone)]
pub struct Dual {
    pub f: f64,
    pub df: f64,
}

impl Dual {
    pub(crate) fn tanh(p0: Dual) -> Dual {
        Dual{f: f64::tanh(p0.f), df: f64::powi(f64::cosh(p0.f),-2)*p0.df, }
    }
}

impl Dual {
    pub(crate) fn cosh(p0: Dual) -> Dual {
        Dual{f: f64::cosh(p0.f), df: f64::sinh(p0.f)*p0.df, }
    }
}

impl From<f64> for Dual {
    fn from(f: f64) -> Dual {
        Dual { f: f, df: 0.0 }
    }
}

impl PartialEq for Dual {
    fn eq(&self, other: &Dual) -> bool {
        self.f == other.f && self.df == other.df
    }
}

impl Ord for Dual {
    fn cmp(&self, other: &Dual) -> Ordering {
        self.f.partial_cmp(&other.f).unwrap()
    }
}
impl Eq for Dual {
    fn assert_receiver_is_total_eq(&self) {
        assert!(self.f.is_nan());
    }
}
impl PartialOrd for Dual {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.f.partial_cmp(&other.f)
    }

    fn lt(&self, other: &Self) -> bool {
        self.f < other.f
    }

    fn le(&self, other: &Self) -> bool {
        self.f <= other.f
    }

    fn gt(&self, other: &Self) -> bool {
        self.f > other.f
    }

    fn ge(&self, other: &Self) -> bool {
        self.f >= other.f
    }
}
impl Dual {
    pub const EPSILON: Dual = Dual { f: 0.0, df: 0.0 };
    pub(crate) fn powi(p0: Dual, p1: i32) -> Dual {
        Dual{
            f: f64::powi(p0.f, p1),
            df: (p1 as f64)*f64::powi(p0.f, p1-1)*p0.df
        }
    }
}

impl Dual {
    pub(crate) fn sqrt(p0: Dual) -> Dual {
        Dual {f:f64::sqrt(p0.f), df:0.5/f64::sqrt(p0.f)*p0.df}
    }
}

impl Dual {
    pub(crate) fn diff(p0: Box<dyn Fn(Dual) -> Dual>, p1:f64) -> f64 {
        let epsilon = Dual::new(p1, 1.0);
        let x = p0(epsilon).df;
        x
    }
}

impl Dual {
    pub fn new(f: f64, df: f64) -> Dual {
        Dual { f, df }
    }
}

impl Add for Dual {
    type Output = Dual;
    fn add(self, rhs:Self) -> Dual {
        Dual{f: self.f+rhs.f, df: self.df+rhs.df}
    }
}

impl Sub for Dual {
    type Output = Dual;
    fn sub(self, rhs:Self) -> Dual {
        Dual{f: self.f-rhs.f, df: self.df-rhs.df}
    }
}

impl Mul for Dual {
    type Output = Dual;
    fn mul(self, rhs:Self) -> Dual {
        Dual{f: self.f*rhs.f, df: self.f*rhs.df+self.df*rhs.f}
    }
}

impl Div for Dual {
    type Output = Dual;
    fn div(self, rhs:Self) -> Dual {
        Dual{f: self.f/rhs.f, df: self.df/rhs.f-self.f*rhs.df/(rhs.f*rhs.f)}
    }
}

impl Add<Dual> for f64 {
    type Output = Dual;

    fn add(self, rhs: Dual) -> Self::Output {
        Dual{f:self+rhs.f, df:rhs.df}
    }
}