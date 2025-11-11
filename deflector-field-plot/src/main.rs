#![allow(refining_impl_trait)]

mod csv_dump;

use nalgebra::Vector4;
use std::env;
use std::str::FromStr;
use deflector_core::warp_drive::{WarpDrive, WarpDriveVBase};
use deflector_core::wd_natario::WarpDriveNatario;
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

fn dump_vx(
    steps: i32,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    file_name: &str,
    wd: &dyn WarpDrive
) {
    dumper(
        steps,
        x_min,
        x_max,
        y_min,
        y_max,
        file_name,
        |v| wd.vx(v)
    );
}

fn dump_v_base(
    steps: i32,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    file_name: &str,
    wd: &dyn WarpDriveVBase
) {
    dumper(
        steps,
        x_min,
        x_max,
        y_min,
        y_max,
        file_name,
        |v| wd.v_base(v)
    );
}

fn dumper<F: FnMut(&Vector4<f64>) -> f64>(
    steps: i32,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    file_name: &str,
    mut function: F
) {
    let x_step = (x_max - x_min) / steps as f64;
    let y_step = (y_max - y_min) / steps as f64;
    let mut dump = CsvDump::new(file_name);

    for i in 0..=steps {
        let x = x_min + x_step * i as f64;
        for j in 0..=steps {
            let y = y_min + y_step * j as f64;
            let val = function(&Vector4::new(0., x, y, 0.));
            dump.add_record(DataPoint3d { x, y, z: 0.0, val });
        }
    }

    dump.dump().expect("Failed to dump CSV");
}

fn print_usage(args: &[String]) {
    eprintln!("Usage: {} <kind> <radius> <sigma> <u> [u0] [k0] [deflector_sigma_pushout] [deflector_sigma_factor] [deflector_back]", args[0]);
}

fn print_usage_ours(args: &[String]) {
    eprintln!("Usage: {} ours <radius> <sigma> <u> <u0> <k0> <deflector_sigma_pushout> <deflector_sigma_factor> <deflector_back>", args[0]);
}

fn print_usage_natario(args: &[String]) {
    eprintln!("Usage: {} natario <radius> <sigma> <u>", args[0]);
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let kind = if args.len() >= 2 {
        args[1].as_str().to_lowercase()
    } else {
        print_usage(&args);
        return;
    };

    let wd: Box<dyn WarpDrive> = match kind.as_str() {
        "ours" if args.len() != 10 => {
            print_usage_ours(&args);
            return;
        },
        "natario" if args.len() != 5 => {
            print_usage_natario(&args);
            return;
        }
        "ours" => {
            let radius = f64::from_str(&args[2]).expect("Failed to parse radius");
            let sigma = f64::from_str(&args[3]).expect("Failed to parse sigma");
            let u = f64::from_str(&args[4]).expect("Failed to parse u");
            let u0 = f64::from_str(&args[5]).expect("Failed to parse u0");
            let k0 = f64::from_str(&args[6]).expect("Failed to parse k0");
            let deflector_sigma_pushout = f64::from_str(&args[7]).expect("Failed to parse deflector_sigma_pushout");
            let deflector_sigma_factor = f64::from_str(&args[8]).expect("Failed to parse deflector_sigma_factor");
            let deflector_back = f64::from_str(&args[9]).expect("Failed to parse deflector_back");

            println!("radius: {}", radius);
            println!("sigma: {}", sigma);
            println!("u: {}", u);
            println!("u0: {}", u0);
            println!("k0: {}", k0);
            println!("deflector_sigma_pushout: {}", deflector_sigma_pushout);
            println!("deflector_sigma_factor: {}", deflector_sigma_factor);
            println!("deflector_back: {}", deflector_back);

            Box::new(WarpDriveOurs::new(radius, sigma, u, u0, k0, deflector_sigma_pushout, deflector_sigma_factor, deflector_back))
        },
        "natario" => {
            let radius = f64::from_str(&args[2]).expect("Failed to parse radius");
            let sigma = f64::from_str(&args[3]).expect("Failed to parse sigma");
            let u = f64::from_str(&args[4]).expect("Failed to parse u");

            println!("radius: {}", radius);
            println!("sigma: {}", sigma);
            println!("u: {}", u);

            Box::new(WarpDriveNatario::new(radius, sigma, u))
        }
        _ => {
            eprintln!("Unknown kind: {}", kind);
            return;
        }
    };

    let x_min = -12.;
    let x_max = 12.;
    let y_min = -12.;
    let y_max = 12.;
    let steps = 24 * 10;

    dump_vx(steps, x_min, x_max, y_min, y_max, "plot_vx.csv", &*wd);

    if let Some(wd) = wd.as_warp_drive_v_base() {
        dump_v_base(steps, x_min, x_max, y_min, y_max, "plot_base.csv", wd);
    }
}
