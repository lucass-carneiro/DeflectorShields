pub fn poly_trans_5(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 {
        y0
    } else if x > (x0 + dx) {
        0.0
    } else {
        ((f64::powi(dx, 2) + 3.0 * dx * (x - x0) + 6.0 * f64::powi(x - x0, 2))
            * f64::powi(dx - x + x0, 3)
            * y0)
            / f64::powi(dx, 5)
    }
}

pub fn d_poly_trans_5_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 || x > (x0 + dx) {
        0.0
    } else {
        (-30.0 * f64::powi(x - x0, 2) * f64::powi(dx - x + x0, 2) * y0) / f64::powi(dx, 5)
    }
}

pub fn poly_trans_7(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 {
        y0
    } else if x > (x0 + dx) {
        0.0
    } else {
        ((f64::powi(dx, 3)
            + 4.0 * f64::powi(dx, 2) * (x - x0)
            + 10.0 * dx * f64::powi(x - x0, 2)
            + 20.0 * f64::powi(x - x0, 3))
            * f64::powi(dx - x + x0, 4)
            * y0)
            / f64::powi(dx, 7)
    }
}

pub fn d_poly_trans_7_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 || x > (x0 + dx) {
        0.0
    } else {
        (-140.0 * f64::powi(x - x0, 3) * f64::powi(dx - x + x0, 3) * y0) / f64::powi(dx, 7)
    }
}

pub fn poly_trans_9(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 {
        y0
    } else if x > (x0 + dx) {
        0.0
    } else {
        ((f64::powi(dx, 4)
            + 5.0 * f64::powi(dx, 3) * (x - x0)
            + 15.0 * f64::powi(dx, 2) * f64::powi(x - x0, 2)
            + 35.0 * dx * f64::powi(x - x0, 3)
            + 70.0 * f64::powi(x - x0, 4))
            * f64::powi(dx - x + x0, 5)
            * y0)
            / f64::powi(dx, 9)
    }
}

pub fn d_poly_trans_9_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    if x < x0 || x > (x0 + dx) {
        0.0
    } else {
        (-630.0 * f64::powi(x - x0, 4) * f64::powi(dx - x + x0, 4) * y0) / f64::powi(dx, 9)
    }
}

const SMALLNESS_TOLERANCE: f64 = 1.0e-2;

pub fn cinf(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    let x0mx = x0 - x;

    // Close to inner radius limit
    if f64::abs(x0mx) < SMALLNESS_TOLERANCE {
        y0
    } else if x0 < x && x < (x0 + dx) {
        let x1 = 1.0 / x0mx;
        let x2 = 1.0 / (x0mx + dx);
        let arg = dx * (x1 + x2);

        let sech_arg = 1.0 / f64::cosh(arg);
        let tanh_arg = f64::tanh(arg);

        (y0 * sech_arg) / (1.0 + tanh_arg + sech_arg)
    } else if x > x0 {
        0.0
    } else {
        y0
    }
}

pub fn d_cinf_dx(x: f64, y0: f64, x0: f64, dx: f64) -> f64 {
    let x0mx = x0 - x;

    // Close to inner radius limit
    if f64::abs(x0mx) < SMALLNESS_TOLERANCE {
        0.0
    } else if x0 < x && x < (x0 + dx) {
        let x0mxpdx = x0mx + dx;

        let x1 = 1.0 / x0mx;
        let x2 = 1.0 / x0mxpdx;
        let arg = (dx / 2.0) * (x1 + x2);

        let sech_arg = 1.0 / f64::cosh(arg);

        -(y0 * dx * (2.0 * x0mx * x0mx + 2.0 * x0mx * dx + dx * dx) * sech_arg * sech_arg)
            / (4.0 * x0mx * x0mx * x0mxpdx * x0mxpdx)
    } else {
        0.0
    }
}
