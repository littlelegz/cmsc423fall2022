use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;
use std::cmp::Ordering;
use serde::{Serialize, Deserialize};
use bincode;

// Struct to serialize
#[derive(Serialize, Deserialize, Debug)]
struct FMindex {
    genome: String,
    bwt: Vec<String>,
    sa: Vec<usize>,
    first_column: Vec<u32>,
    tally: Vec<Vec<u32>>,
}

fn main() {
    // Read input args
    let args: Vec<_> = env::args().collect();
    if args.len() != 3 {
        println!("Invalid parameters");
        exit(1);
    }
    let reference: &str = &args[1];
    let output: &str = &args[2];
    // Reading genome
    let reads: &str = &read_gen(reference);
    //let reads: &str = &reads_forward.chars().rev().collect::<String>();
    // Define suffix saay
    let mut sa: Vec<usize> = Vec::new();
    // Fill with suffix indexes 
    for i in 0..reads.len() {
        sa.push(i);
    }

    // Quicksort suffix sa indices
    quick_sort::<usize>(&mut sa, reads);
    //println!("Sorted sa: {:?}", sa);

    /* Converting sa to BWT 
       Also counting the tally in form <$: ?, A: ?, C: ?, G: ?, T: ?> and saving the first column */
    let mut bwt: Vec<String> = Vec::new();
    let mut tally: Vec<Vec<u32>> = Vec::new();
    let mut tally_counter: Vec<u32> = vec![0,0,0,0,0];

    for index in &sa {
        let mut curr: String = "".to_string();
        if *index == 0 as usize {
            bwt.push("$".to_string());
        } else {
            curr = reads[index - 1 as usize..*index].to_string();
            bwt.push(reads[index - 1 as usize..*index].to_string());
        }

        // Tallying occ of bases in BWT
        match curr.as_str() {
            "A" => tally_counter[1] += 1,
            "C" => tally_counter[2] += 1,
            "G" => tally_counter[3] += 1,
            "T" => tally_counter[4] += 1,
            &_ => tally_counter[0] += 1,
        }
        tally.push(tally_counter.clone())
    }

    //println!("BWT: {:?}", bwt);
    println!("Tally: {:?}", tally_counter);

    // Creating output file
    println!("[*] Creating output file: {}", output);
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 

    // Serialize data
    let out_struct = FMindex {
        genome: reads.to_string(),
        bwt: bwt,
        sa: sa,
        first_column: tally_counter,
        tally: tally,
    };

    let serial_data: Vec<u8> = bincode::serialize(&out_struct).unwrap();

    // Write to file
    println!("[*] Writing data to binary file");
    write_output.write_all(&serial_data);
}

// Maybe implement custom comparator for efficiency?
pub fn compare (curr: &str, comp: &str) -> Ordering {
    return curr.cmp(comp);
}

// Quicksort recursion
pub fn quick_sort<T: Ord>(sa: &mut Vec<usize>, reads: &str) {
    let len = sa.len();
    _quick_sort::<T>(sa, 0, (len - 1) as isize, reads);
}

fn _quick_sort<T: Ord>(sa: &mut Vec<usize>, low: isize, high: isize, reads: &str) {
    if low < high {
        let p = partition::<T>(sa, low, high, reads);
        _quick_sort::<T>(sa, low, p - 1, reads);
        _quick_sort::<T>(sa, p + 1, high, reads);
    }
}

fn partition<T: Ord>(sa: &mut Vec<usize>, low: isize, high: isize, reads: &str) -> isize {
    let pivot = high as usize;
    let mut store_index = low - 1;
    let mut last_index = high;

    loop {
        store_index += 1;
        while compare(&reads[sa[store_index as usize]..], &reads[sa[pivot]..]) == Ordering::Less  {
            store_index += 1;
        }
        last_index -= 1;
        while last_index >= 0 && compare(&reads[sa[last_index as usize]..], &reads[sa[pivot]..]) == Ordering::Greater {
            last_index -= 1;
        }
        if store_index >= last_index {
            break;
        } else {
            sa.swap(store_index as usize, last_index as usize);
        }
    }
    sa.swap(store_index as usize, pivot as usize);
    store_index
}

pub fn read_gen(path: &str) -> String {
    let genome_file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(genome_file);
    let mut gen: String = "".to_owned();
    
    for line in reader.lines() {
        if let Ok(curr) = line {
            match curr.chars().next() {
                Some(firstchar) => {
                    if firstchar == '>' {
                        continue;
                    } else {
                        gen.push_str(&curr);
                    }
                }
                None => {
                    break;
                }
            }
            
        }
    }
    gen.push_str("$");
    println!("[*] Read genome into gen, appended sentinel, size is: {}", gen.chars().count());

    return gen;
}
