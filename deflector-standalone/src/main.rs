use std::env;

use deflector_core::{evolve, multi_ids, output, params, warp_drive::WarpDrive};

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

fn integrate<Drive: WarpDrive>(output_file_name: &str, param_file_name: &str, warp_drive: &Drive) {
    /*
     * We have to read the parameter file again for borrow checker reasons. This could potentially
     * be avoided by implementing Clone for Params, but that cascades into a lot of changes that
     * I don't want to make rn to not break the current implementation of the simulator
     */
    let par = params::read_params(param_file_name).unwrap();

    // Initialize particle state vectors.
    let mut states =
        multi_ids::make_initial_data(par.initial_data, &par.normalize_as, warp_drive).unwrap();

    // Compute the number of time steps
    let num_time_steps =
        (par.time_integration.final_time / par.time_integration.time_step) as usize;

    for i in 0..=num_time_steps {
        let t = (i as f64) * par.time_integration.time_step;
        log::info!("Integrating step {}/{}, t = {}", i, num_time_steps, t);

        // Do output
        if (i % par.time_integration.out_every.unwrap_or(1usize)) == 0 {
            let mut out_file = output::IpcMultiFile::new();

            // Loop over particles
            for particle_idx in 0..states.len() {
                out_file.append(particle_idx, i, t, &states[particle_idx]);

                evolve::rk4_step(
                    t,
                    par.time_integration.time_step,
                    warp_drive,
                    &mut states[particle_idx],
                );
            }

            out_file.write(output_file_name, i);
        } else {
            // Loop over particles
            for particle_idx in 0..states.len() {
                evolve::rk4_step(
                    t,
                    par.time_integration.time_step,
                    warp_drive,
                    &mut states[particle_idx],
                );
            }
        }
    }
}

fn main() {
    init_logger();

    // Parse cmd line
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        println!("Usage: {} <parameter-file> <output-file>", args[0]);
        return;
    }

    let param_file_name = &args[1];
    let output_file_name = &args[2];

    // Deserialize parameter file and
    let par = params::read_params(param_file_name).unwrap();

    // Integrate warp drive based on the type of drive specified in the parameter file
    match par.warp_drive {
        params::WarpDriveType::Natario(warp_drive_natario) => {
            integrate(output_file_name, param_file_name, &warp_drive_natario)
        }
        params::WarpDriveType::Ours(warp_drive_ours) => {
            integrate(output_file_name, param_file_name, &warp_drive_ours)
        }
    };
}
