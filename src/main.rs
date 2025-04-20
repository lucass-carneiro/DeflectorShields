use metrics::Metric;

mod errors;
mod evolve;
mod logging;
mod metrics;
mod minkowski;
mod params;

fn main() {
    logging::init_logger();

    // Parse args
    let cmd_args: Vec<String> = std::env::args().collect();

    if cmd_args.len() != 2 {
        log::error!("Usage: {} <parameter-file>", cmd_args[0]);
        return;
    }

    // Parse params
    log::info!("Reading parameter file \"{}\"", cmd_args[1]);

    let par = match params::read_params(&cmd_args[1]) {
        Ok(o) => o,
        Err(e) => {
            log::error!("{}", e.to_string());
            return;
        }
    };

    // Time step data
    let tf = par.time_integration.tf;
    let dt = par.time_integration.dt;
    let nt = (f64::ceil(tf / dt) as u64) + 1;

    // Normalize and evolve
    match par.spacetime {
        params::Spacetime::Minkowski => {
            let metric = minkowski::Minkowski {};

            log::info!("Normalizing initial p_0");
            let mut state = evolve::normalize(par, &metric).unwrap();

            log::info!("Normalized initial p_0 is {}", state.p[0]);

            for i in 0..nt {
                let t = (i as f64) * dt;
                log::info!(
                    "({t}, {}, {}, {}, {})",
                    state.x[0],
                    state.x[1],
                    state.x[2],
                    state.x[3]
                );
                evolve::rk4_step(&mut state, &metric, dt);
            }
        }

        params::Spacetime::Alcubierre => {
            // TODO
            return;
        }
    }
}
