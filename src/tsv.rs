use std::fs::File;
use std::io::Write;

use crate::evolve::StateVector;

pub fn make_file(file_name: &str) -> File {
    let mut file = File::create(file_name).unwrap();
    writeln!(
        &mut file,
        "# 1: Iteration 2:Lambda 3:t 4:x 5:y 6:z 7:pt 8:px 9:py 10:pz"
    )
    .unwrap();
    file
}

pub fn append_data(file: &mut File, iteration: u64, lambda: f64, state_vector: &StateVector) {
    writeln!(
        file,
        "{}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}\t{:.16}",
        iteration,
        lambda,
        state_vector.x[0],
        state_vector.x[1],
        state_vector.x[2],
        state_vector.x[3],
        state_vector.p[0],
        state_vector.p[1],
        state_vector.p[2],
        state_vector.p[3],
    )
    .unwrap();
}
