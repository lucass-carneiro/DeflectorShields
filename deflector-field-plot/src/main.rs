#![allow(refining_impl_trait)]

mod csv_dump;

use nalgebra::Vector4;
use std::env;
use std::str::FromStr;
use deflector_core::warp_drive::WarpDrive;
use deflector_core::wd_ours::WarpDriveOurs;
use crate::csv_dump::{CsvDump, CsvRecord};

struct DataPoint3d {
    x: f64,
    y: f64,
    z: f64,
    val: f64
}

impl CsvRecord<4> for DataPoint3d {
    fn get_csv_header() -> [&'static str; 4] {
        ["x", "y", "z", "val"]
    }

    fn get_csv_record(&self) -> [f64; 4] {
        [self.x, self.y, self.z, self.val]
    }
}

fn dumper(
    steps: i32,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    fname: &str,
    base: bool,
    wd: &WarpDriveOurs
) {
    let x_step = (x_max - x_min) / steps as f64;
    let y_step = (y_max - y_min) / steps as f64;
    let mut dump = CsvDump::new(fname);

    for i in 0..=steps {
        let x = x_min + x_step * i as f64;
        for j in 0..=steps {
            let y = y_min + y_step * j as f64;

            if base {
                let val = wd.v_base(&Vector4::new(0., x, y, 0.));
                dump.add_record(DataPoint3d { x, y, z: 0.0, val });
            } else {
                let val = wd.vx(&Vector4::new(0., x, y, 0.));
                dump.add_record(DataPoint3d { x, y, z: 0.0, val });
            }
        }
    }

    dump.dump().expect("Failed to dump CSV");
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let radius = f64::from_str(&args[1]).unwrap();
    let sigma = f64::from_str(&args[2]).unwrap();
    let u = f64::from_str(&args[3]).unwrap();
    let u0 = f64::from_str(&args[4]).unwrap();
    let k0 = f64::from_str(&args[5]).unwrap();
    let deflector_sigma_pushout = f64::from_str(&args[6]).unwrap();
    let deflector_sigma_factor = f64::from_str(&args[7]).unwrap();
    let deflector_back = f64::from_str(&args[8]).unwrap();
    println!("radius: {}", radius);
    println!("sigma: {}", sigma);
    println!("u: {}", u);
    println!("u0: {}", u0);
    println!("k0: {}", k0);
    println!("deflector_sigma_pushout: {}", deflector_sigma_pushout);
    println!("deflector_sigma_factor: {}", deflector_sigma_factor);
    println!("deflector_back: {}", deflector_back);
    let wd = WarpDriveOurs::new(radius, sigma, u, u0, k0, deflector_sigma_pushout, deflector_sigma_factor, deflector_back);

    let x_min = -12.;
    let x_max = 12.;
    let y_min = -12.;
    let y_max = 12.;

    let steps = 24 * 10;

    dumper(steps, x_min, x_max, y_min, y_max, "plot_base.csv", true, &wd);
    dumper(steps, x_min, x_max, y_min, y_max, "plot_vx.csv", false, &wd);

    // {
    //     let mut dump = CsvDump::new("plot_base.csv");
    //
    //     for i in 0..=steps {
    //         let x = x_min + x_step * i as f64;
    //         for j in 0..=steps {
    //             let y = y_min + y_step * j as f64;
    //
    //             let val = wd.v_base(&Vector4::new(0., x, y, 0.));
    //
    //             dump.add_record(DataPoint3d { x, y, z: 0.0, val });
    //         }
    //     }
    //
    //     dump.dump().expect("Failed to dump CSV");
    // }
}
