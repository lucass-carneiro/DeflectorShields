#![allow(refining_impl_trait)]

mod csv_dump;

use nalgebra::Vector4;
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

fn main() {
    let wd = WarpDriveOurs::new(
        4.,
        4.,
        0.5,
        0.5,
        0.1,
        1.0,
        1.1,
        0.5,
    );

    let x_min = -12.;
    let x_max = 12.;
    let y_min = -12.;
    let y_max = 12.;

    let steps = 24 * 10;
    let x_step = (x_max - x_min) / steps as f64;
    let y_step = (y_max - y_min) / steps as f64;

    let mut dump = CsvDump::new("plot.csv");

    for i in 0..=steps {
        let x = x_min + x_step * i as f64;
        for j in 0..=steps {
            let y = y_min + y_step * j as f64;

            let val = wd.vx(&Vector4::new(0., x, y, 0.));

            dump.add_record(DataPoint3d {
                x,
                y,
                z: 0.0,
                val
            });
        }
    }

    dump.dump().expect("Failed to dump CSV");
}
