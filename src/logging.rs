pub fn init_logger() {
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
