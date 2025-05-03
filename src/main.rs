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

fn single_particle_mode(output_file_name: &str, par: params::Params) {
    log::info!("Single particle mode");
    let id = par.single_particle_id.unwrap();

    let warp_drive_solution: Box<dyn WarpDriveHamiltonian> = match par.warp_drive_solution {
        params::WarpDriveSolution::Alcubierre(sol) => Box::new(sol),
        params::WarpDriveSolution::AlcubierreSharp(sol) => Box::new(sol),
    };

    let mut state = warp_drive_solution
        .make_normalized_state(
            0.0,
            id.x0,
            id.y0,
            id.z0,
            id.px0,
            id.py0,
            id.pz0,
            &par.normalize_as,
        )
        .unwrap();

    let mut out_file = output::IpcFile::new();

    let nlambda = (par.affine_data.lambda_max / par.affine_data.dlambda) as u64;

    for i in 0..=nlambda {
        let t = (i as f64) * par.affine_data.dlambda;
        log::info!("Integrating step {}/{}, t = {}", i, nlambda, t);

        out_file.append(i, t, &state);

        evolve::rk4_step(par.affine_data.dlambda, &warp_drive_solution, &mut state);
    }

    out_file.write(output_file_name);

    log::info!("Single particle mode finished");
}

fn multi_particle_mode(output_file_name: &str, par: params::Params) {
    log::info!("Multiple particles mode");

    let id = par.multi_particle_id.unwrap();

    let warp_drive_solution: Box<dyn WarpDriveHamiltonian> = match par.warp_drive_solution {
        params::WarpDriveSolution::Alcubierre(sol) => Box::new(sol),
        params::WarpDriveSolution::AlcubierreSharp(sol) => Box::new(sol),
    };

    let mut states = multi_ids::make_multi_id(id, &par.normalize_as, &warp_drive_solution).unwrap();

    let nlambda = (par.affine_data.lambda_max / par.affine_data.dlambda) as usize;

    for i in 0..=nlambda {
        let t = (i as f64) * par.affine_data.dlambda;
        log::info!("Integrating step {}/{}, t = {}", i, nlambda, t);

        if (i % par.affine_data.out_every.unwrap_or(1usize)) == 0 {
            let mut out_file = output::IpcMultiFile::new();

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

    if par.single_particle_id.is_some() && par.multi_particle_id.is_none() {
        single_particle_mode(output_file_name, par);
    } else if par.multi_particle_id.is_some() && par.single_particle_id.is_none() {
        multi_particle_mode(output_file_name, par);
    } else {
        log::error!(
            "Cannot operate when both multi and single particle initial conditions are provided"
        );
    }
}
