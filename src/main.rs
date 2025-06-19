use std::env;

mod errors;
mod evolve;
mod multi_ids;
mod output;
mod params;
mod types;
mod warp_drive_hamiltonian;

mod warp_drive_alcubierre;
mod warp_drive_alcubierre_sharp;

use warp_drive_hamiltonian::WarpDriveHamiltonian;

fn init_logger() {
    use env_logger::{Builder, Env};
    use std::io::Write;

    let env = Env::default().filter_or("DS_LOG_LEVEL", "info");

    Builder::from_env(env)
        .format(|buf, record| {
            writeln!(
                buf,
                "[Deflector Shield {} {}]: {}",
                buf.timestamp(),
                record.level(),
                record.args()
            )
        })
        .init();
}

fn main() {
    init_logger();

    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        println!("Usage: {} <parameter-file> <output-file>", args[0]);
        return;
    }

    let param_file_name = &args[1];
    let output_file_name = &args[2];

    let par = params::read_params(param_file_name).unwrap();

    // Get warp drive solution to use
    let warp_drive_solution: Box<dyn WarpDriveHamiltonian> = match par.warp_drive_solution {
        params::WarpDriveSolution::Alcubierre(sol) => Box::new(sol),
        params::WarpDriveSolution::AlcubierreSharp(sol) => Box::new(sol),
    };

    // Initialize particle state vectors. The ship is allways the last particle
    let mut states = multi_ids::make_multi_id(
        par.multi_particle_id,
        &par.normalize_as,
        par.ship_speed,
        &warp_drive_solution,
    )
    .unwrap();

    // Compute the number of time steps
    let nlambda = (par.affine_data.lambda_max / par.affine_data.dlambda) as usize;

    // Time stepping loop
    for i in 0..=nlambda {
        let t = (i as f64) * par.affine_data.dlambda;
        log::info!("Integrating step {}/{}, t = {}", i, nlambda, t);

        // Do output
        if (i % par.affine_data.out_every.unwrap_or(1usize)) == 0 {
            let mut out_file = output::IpcMultiFile::new();

            // Loop over particles
            for particle_idx in 0..states.len() {
                out_file.append(particle_idx, i, t, &states[particle_idx]);

                evolve::rk4_step(
                    par.affine_data.dlambda,
                    &warp_drive_solution,
                    &mut states[particle_idx],
                );
            }

            out_file.write(&output_file_name, i);
        } else {
            // Loop over particles
            for particle_idx in 0..states.len() {
                evolve::rk4_step(
                    par.affine_data.dlambda,
                    &warp_drive_solution,
                    &mut states[particle_idx],
                );
            }
        }
    }
}
