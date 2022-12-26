use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;
use std::cmp::{Ordering, min};

// Struct to serialize
#[derive(Serialize, Deserialize, Debug)]
pub struct GenSA {
    genome: String,
    sa: Vec<usize>
}

#[derive(Debug, Clone)]
pub struct SaQ {
    header: String,
    query: String,
    hits: Vec<usize>,
    char_cmp_lb: u32,
    char_cmp_ub: u32
}

fn main () {
    let args: Vec<_> = env::args().collect();
    if args.len() != 5 {
        println!("Invalid parameters");
        exit(1);
    }
    let index: &str = &args[1];
    let queries: &str = &args[2];
    let query_mode: &str = &args[3];
    let output: &str = &args[4];
    // Reading sa
    let sa_struct: GenSA = read_sa_struct(index);

    // Reading queries
    let mut query_arr: Vec<SaQ> = read_queries(queries);
    println!("[*] Read queries, count: {}", query_arr.len());

    // Fulfill query
    if query_mode == "naive" {
        for (idx, query) in query_arr.clone().iter().enumerate() {
            let tuple = binary_search_wrapper(&query.query, sa_struct.sa.clone(), &sa_struct.genome);
            query_arr[idx].hits = tuple.0;
            query_arr[idx].char_cmp_lb = tuple.1;
            query_arr[idx].char_cmp_ub = tuple.2;
            //println!("{:?}", query);
        }
    } else if query_mode == "simpaccel" { //Mode is simpaccel
        for (idx, query) in query_arr.clone().iter().enumerate() {
            let tuple = simple_accel_wrapper(&query.query, sa_struct.sa.clone(), &sa_struct.genome);
            query_arr[idx].hits = tuple.0;
            query_arr[idx].char_cmp_lb = tuple.1;
            query_arr[idx].char_cmp_ub = tuple.2;
            //println!("{:?}", query);
        }
    } else {
        println!("Input query mode of: {} was invalid", query_mode);
        return 
    }

    // Output the query struct to output file
    // println!("{:?}", query_arr);

    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    for query in query_arr {
        if let Err(e) = writeln!(write_output, "{}\t{}\t{}\t{}\t{}", query.header, query.char_cmp_lb, query.char_cmp_ub, query.hits.len(),
                                 query.hits.iter().map( |&id| id.to_string() + "\t").collect::<String>()) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}

pub fn simple_accel_wrapper(query: &str, sa: Vec<usize>, genome: &str) -> (Vec<usize>, u32, u32) {
    let query_low: &str = &format!("{}{}",query,"#");
    let query_up: &str = &format!("{}{}",query,"}");

    // Perform lower search
    let lower = simple_accel(query_low, sa.clone(), genome);
    let _lower_string: &str = &genome[sa[lower.0 as usize]..];
    
    // Perform upper search
    let upper = simple_accel(query_up, sa.clone(), genome);
    let _upper_string: &str = &genome[sa[upper.0 as usize]..];
    println!("[*] Found lower at: {}, Upper at: {}, Sa lower: {:?}, Sa Upper: {:?}", lower.0, upper.0, sa[lower.0 as usize], sa[upper.0 as usize]);
    
    return (sa[lower.0 as usize..upper.0 as usize].to_vec(), lower.1, upper.1);
}

pub fn simple_accel(target: &str, sa: Vec<usize>, genome: &str) -> (u32, u32) {
    let mut low: i32 = 0;
    let mut high: i32 = sa.len() as i32;
    let mut comparisons: u32 = 0;
    let mut result_low: u32 = 0;
    let mut result_high: u32 = 0;
 
    loop {
        let middle: i32 = (high + low) / 2;
        let current = &genome[sa[middle as usize]..];
        let result: (Ordering, u32) = compare(&target[min(result_low, result_high) as usize..], &current[min(result_low, result_high) as usize..]);
        comparisons += result.1;

        if result.0 == Ordering::Less {
            if middle == low + 1 {
                return (middle as u32, comparisons);
            }
            high = middle;
            result_high = (result.1 - 1) + min(result_low, result_high);
        }
        if result.0 == Ordering::Greater {
            if middle == high - 1 {
                return (high as u32, comparisons);
            }
            low = middle;
            result_low = (result.1 - 1) + min(result_low, result_high);
        }
    }
}

// Binary search into sa
pub fn binary_search_wrapper(query: &str, sa: Vec<usize>, genome: &str) -> (Vec<usize>, u32, u32) {
    let query_low: &str = &format!("{}{}",query,"#");
    let query_up: &str = &format!("{}{}",query,"}");

    // Perform lower search
    let lower = binary_search(query_low, sa.clone(), genome);
    let _lower_string: &str = &genome[sa[lower.0 as usize]..];
    
    // Perform upper search
    let upper = binary_search(query_up, sa.clone(), genome);
    let _upper_string: &str = &genome[sa[upper.0 as usize]..];
    println!("[*] Found lower at: {}, Upper at: {}, Sa lower: {:?}, Sa Upper: {:?}", lower.0, upper.0, sa[lower.0 as usize], sa[upper.0 as usize]);
    
    return (sa[lower.0 as usize..upper.0 as usize].to_vec(), lower.1, upper.1);
}

pub fn binary_search(target: &str, sa: Vec<usize>, genome: &str) -> (u32, u32) {
    let mut low: i32 = 0;
    let mut high: i32 = sa.len() as i32;
    let mut comparisons: u32 = 0;
 
    loop {
        let middle: i32 = (high + low) / 2;
        let current = &genome[sa[middle as usize]..];
        let result = compare(target, &current);
        comparisons += result.1;

        if result.0 == Ordering::Less {
            if middle == low + 1 {
                return (middle as u32, comparisons);
            }
            high = middle
        }
        if result.0 == Ordering::Greater {
            if middle == high - 1 {
                return (high as u32, comparisons);
            }
            low = middle
        }
    }
}

// Maybe implement custom comparator
pub fn compare(curr: &str, comp: &str) -> (Ordering, u32) {
    let mut order: Ordering = Ordering::Equal;
    let iter = curr.chars().zip(comp.chars());
    let mut comp_count: u32 = 0;
    
    for (a, b) in iter {
        comp_count += 1;
        if a != b {
            order = a.cmp(&b);
            break;
        }
    }

    //println!("Compared {} with {}, comp_count was: {}", curr, &comp[0..curr.len()], comp_count);
    return (order, comp_count);
}

// Modify to return vec<str> of queries
pub fn read_queries(path: &str) -> Vec<SaQ> {
    let query_file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(query_file);
    let mut query_arr: Vec<SaQ> = Vec::new();
    let mut query: String = "".to_owned();
    let mut header: String = "".to_owned();
    
    for line in reader.lines() {
        if let Ok(mut curr) = line {
            match curr.chars().next() {
                Some(firstchar) => {
                    if firstchar == '>' {
                        curr.remove(0);
                        if query.is_empty() {
                            header = curr.to_string();
                            continue;
                        }
                        query_arr.push(
                            SaQ {
                                header: header.to_string(),
                                query: query.to_string(),
                                hits: Vec::new(),
                                char_cmp_lb: 0,
                                char_cmp_ub: 0
                            }
                        );
                        header = curr.to_string();
                        query = "".to_string();
                    } else {
                        query.push_str(&curr);
                    }
                }
                None => {
                    break;
                }
            }
            
        }
    }

    // Push last element
    query_arr.push({
        SaQ {
            header: header.to_string(),
            query: query.to_string(),
            hits: Vec::new(),
            char_cmp_lb: 0,
            char_cmp_ub: 0
        }
    });
    return query_arr;
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
