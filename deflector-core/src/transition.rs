pub fn poly_trans(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 {
        y0
    } else if x > (x0 + dx) {
        0.0
    } else {
        let xmx0 = x - x0;
        ((dx * dx + 3.0 * dx * xmx0 + 6.0 * xmx0 * xmx0) * f64::powi(dx - x + x0, 3) * y0)
            / f64::powi(dx, 5)
    }
}

pub fn d_poly_trans_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 || x > (x0 + dx) {
        0.0
    } else {
        let a = x - x0;
        let b = dx - x + x0;
        -(30.0 * a * a * b * b * y0) / f64::powi(dx, 5)
    }
}

pub fn cinf(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x0 < x && x < (x0 + dx) {
        let x1 = 1.0 / (x0 - x);
        let x2 = 1.0 / (x1 + dx);
        let x3 = dx * (x1 + x2);
        y0 / (1.0 + f64::exp(x3))
    } else if x > x0 {
        0.0
    } else {
        y0
    }
}

pub fn d_cinf_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x0 < x && x < (x0 + dx) {
        let x1 = x0 - x;
        let x2 = x1 + dx;
        let x3 = (dx / 2.0) * (1.0 / x1 + 1.0 / x2);
        let x3 = 1.0 / f64::cosh(x3);
        let x4 = x3 * x3;
        -(y0 * dx * (2.0 * x1 * x1 + 2.0 * x1 * dx + dx * dx) * x4) / (4.0 * x1 * x1 * x2 * x2)
    } else {
        0.0
    }
}
