use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::cmp;
use std::io::{prelude::*, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;

#[derive(Serialize, Deserialize, Debug)]
pub struct GenSA {
    genome: String,
    sa: Vec<usize>
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
    // Reading sa struct
    let sa_struct: GenSA = read_sa_struct(input);
    let sa_len: usize = sa_struct.sa.len();
    // Calculating LCP1
    let lcp1: Vec<i32> = calc_lcp1(&sa_struct.sa, &sa_struct.genome);
    let mut lcp1_copy: Vec<i32> = lcp1.clone();
    println!("[*] Calculated lcp1, length {}", lcp1.len());
    let lcp1_max: i32 = *lcp1.iter().max().expect("lcp not empty");
    let lcp1_med: f64 = get_median(&mut lcp1_copy);
    let lcp1_avg: f64 = lcp1.iter().sum::<i32>() as f64 / lcp1.len() as f64;
    println!("[*] Calculated mean: {:5}, median: {}, maximum {}", lcp1_avg, lcp1_med, lcp1_max);

    let mut sample_str: String = "".to_string();
    // Add first entry
    sample_str.push_str(&sa_struct.sa.get(0 as usize).expect("out of range of sa").to_string());
    let mut i: i32 = 1;
    while i*sample_rate < sa_len as i32 {
        sample_str.push_str("\t");
        sample_str.push_str(&sa_struct.sa.get((i*sample_rate) as usize).expect("out of range of sa").to_string());
        i += 1;
    }
    println!("[*] Calculated sample arr of size {}", sample_str.len());
     
    // Creating output file
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    if let Err(e) = writeln!(write_output, "{}\n{}\n{}\n{}", lcp1_avg, lcp1_med, lcp1_max, sample_str) {
        eprintln!("Couldn't write to file: {}", e);
    }

}

pub fn get_median (arr: &mut Vec<i32>) -> f64 {
    arr.sort();
    if arr.len() % 2 == 0 {
        let mid = (arr.len() - 1) / 2;
        let sum: f64 = (arr.get(mid).expect("out of range of arr") + arr.get(mid + 1).expect("out of range of arr")) as f64;
        return sum / 2 as f64;
    } else {
        let mid = (arr.len() - 1) / 2;
        return *arr.get(mid).expect("out of range of arr") as f64;
    }
}

pub fn get_lcp (strings: Vec<&str>) -> i32 {
    let curr = strings[0];
    let currbytes = curr.as_bytes(); // Why as bytes? Why not just .iter()?
    let mut len = curr.len();
    for str in &strings[1..] {
        len = cmp::min(len, // Take min overlap b/t curr and others
                       str.as_bytes()
                          .iter() // Iterate other
                          .zip(currbytes) // Also iterate over curr
                          .take_while(|&(a, b)| a == b) // Only count if tuple exists, and chars are equal
                          .count()); // Counting
    }
    return len as i32;
}

pub fn calc_lcp1 (arr: &Vec<usize>, genome: &str) -> Vec<i32> {
    let mut lcp1 = Vec::new();
    for i in 0..arr.len()-1 {
        lcp1.push(get_lcp(vec![&genome[*arr.get(i).expect("index out of sa")..], &genome[*arr.get(i+1).expect("index out of sa")..]]))
    }
    return lcp1;
}

pub fn read_sa_struct(path: &str) -> GenSA {
    // Reading sa binary file
    let file = File::open(path);
    // read the same file back into a Vec of bytes
    let mut buffer = Vec::<u8>::new();
    file.expect("Unable to open file").read_to_end(&mut buffer);

    let sa_struct: GenSA = bincode::deserialize(&buffer).unwrap();
    println!("[*] Read genome of size: {}, sa of size: {}", sa_struct.genome.len(), sa_struct.sa.len());
    return sa_struct;
}
