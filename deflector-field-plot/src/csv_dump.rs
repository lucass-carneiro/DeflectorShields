use std::fmt::Display;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

pub(super) trait CsvRecord<const FIELDS: usize> {
    fn get_csv_header() -> [&'static str; FIELDS];
    fn get_csv_record(&self) -> [impl Display; FIELDS];

    fn get_csv_header_str() -> String {
        Self::get_csv_header().join(",")
    }

    fn get_csv_record_str(&self) -> String {
        self.get_csv_record().map(|i| i.to_string()).join(",")
    }
}

pub(super) struct CsvDump<R: CsvRecord<FIELDS>, const FIELDS: usize> {
    file_path: PathBuf,
    records: Vec<R>
}

impl<R: CsvRecord<F>, const F: usize> CsvDump<R, F> {
    pub(super) fn new(file_path: impl Into<PathBuf>) -> Self {
        Self {
            file_path: file_path.into(),
            records: Vec::new()
        }
    }

    pub(super) fn add_record(&mut self, record: R) {
        self.records.push(record);
    }

    pub(super) fn dump(&self) -> Result<(), std::io::Error> {
        let mut file =
            OpenOptions::new()
                        .create(true)
                        .write(true)
                        .truncate(true)
                        .open(&self.file_path)?;
        
        let mut writer = BufWriter::new(&mut file);

        writeln!(&mut writer, "{}", R::get_csv_header_str())?;

        for record in &self.records {
            writeln!(&mut writer, "{}", record.get_csv_record_str())?;
        }

        Ok(())
    }
}