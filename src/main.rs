use std::env;

mod alcubierre;
mod errors;
mod evolve;
mod output;
mod params;
mod types;

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

    let par = params::read_params(&args[1]).unwrap();

    match par.single_particle_id {
        Some(id) => {
            log::info!("Single particle mode");

            let mut state = par
                .alcubierre_data
                .make_normalized_state(
                    id.x0,
                    id.y0,
                    id.z0,
                    id.px0,
                    id.py0,
                    id.pz0,
                    par.normalize_as,
                )
                .unwrap();

            let mut out_file = output::TSVFile::new(&args[2]);

            let nlambda = (par.affine_data.lambda_max / par.affine_data.dlambda) as u64;

            for i in 0..=nlambda {
                let t = (i as f64) * par.affine_data.dlambda;
                log::info!("Integrating step {}/{}, t = {}", i, nlambda, t);

                out_file.append(i, t, &state);

                evolve::rk4_step(par.affine_data.dlambda, &par.alcubierre_data, &mut state);
            }

            log::info!("Single particle mode finished");
        }

        None => log::info!("Multi particle mode not implemented yet"),
    }
}
