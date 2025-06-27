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
