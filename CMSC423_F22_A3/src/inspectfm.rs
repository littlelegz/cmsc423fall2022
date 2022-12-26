use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;

#[derive(Serialize, Deserialize, Debug)]
pub struct FMindex {
    bwt: Vec<String>,
    sa: Vec<usize>,
    first_column: Vec<u32>,
    tally: Vec<Vec<u32>>,
}

fn main() {
    // Read input args
    let args: Vec<_> = env::args().collect();
    if args.len() != 4 {
        println!("Invalid parameters");
        exit(1);
    }
    let input: &str = &args[1];
    let sample_rate: i32 = args[2].parse().unwrap();
    let output: &str = &args[3]; 
    // Reading fm struct
    let fm_struct: FMindex = read_fm_struct(input);

    // Sampling FM index
    let mut spot_a: Vec<String> = Vec::new();
    let mut spot_c: Vec<String> = Vec::new();
    let mut spot_g: Vec<String> = Vec::new();
    let mut spot_t: Vec<String> = Vec::new();
    // Add first entry

    //println!("[*] Calculated tally arr of size {}", fm_struct.tally.len());
    let mut i: i32 = 0;
    while i*sample_rate < fm_struct.tally.len() as i32 {
        spot_a.push(fm_struct.tally[(i*sample_rate) as usize][1].to_string());
        spot_c.push(fm_struct.tally[(i*sample_rate) as usize][2].to_string());
        spot_g.push(fm_struct.tally[(i*sample_rate) as usize][3].to_string());
        spot_t.push(fm_struct.tally[(i*sample_rate) as usize][4].to_string());
        i += 1;
    }
    println!("[*] Calculated sample arr of size {}", spot_a.len());
     
    // Creating output file
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    if let Err(e) = writeln!(write_output, "{}\t{}\t{}\t{}\t{}\n{}\n{}\n{}\n{}\n{}", fm_struct.first_column[0],
                             fm_struct.first_column[1], fm_struct.first_column[2], fm_struct.first_column[3],
                             fm_struct.first_column[4], fm_struct.bwt.join(""), spot_a.join("\t"), spot_c.join("\t"),
                             spot_g.join("\t"), spot_t.join("\t")) {
        eprintln!("Couldn't write to file: {}", e);
    }

}

pub fn read_fm_struct(path: &str) -> FMindex {
    // Reading sa binary file
    let file = File::open(path);
    // read the same file back into a Vec of bytes
    let mut buffer = Vec::<u8>::new();
    file.expect("Unable to open file").read_to_end(&mut buffer);

    let fm_struct: FMindex = bincode::deserialize(&buffer).unwrap();
    println!("[*] Read bwt of size: {}, sa of size: {}", fm_struct.bwt.len(), fm_struct.sa.len());
    return fm_struct;
}
