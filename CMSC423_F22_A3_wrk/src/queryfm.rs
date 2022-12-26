use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;
use std::cmp::{Ordering};

// Struct to serialize
#[derive(Serialize, Deserialize, Debug)]
pub struct FMindex {
    bwt: Vec<String>,
    sa: Vec<usize>,
    first_column: Vec<u32>,
    tally: Vec<Vec<u32>>,
} 

#[derive(Debug, Clone)]
pub struct FMQ {
    header: String,
    query: String,
    match_len: u32,
    hits: Vec<usize>
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
    let fm_struct: FMindex = read_fm_struct(index);

    // Reading queries
    let mut query_arr: Vec<FMQ> = read_queries(queries);
    println!("[*] Read queries, count: {}", query_arr.len());

    // Fulfill query
    println!("[*] Fulfilling queries...");
    if query_mode == "complete" {
        for (idx, query) in query_arr.clone().iter().enumerate() {
            let q_out: (Vec<usize>, u32) = complete_search(&query.query, &fm_struct);
            query_arr[idx].hits = q_out.0;
            query_arr[idx].match_len = q_out.1;
        }
    } else if query_mode == "partial" { //Mode is simpaccel
        for (idx, query) in query_arr.clone().iter().enumerate() {
            let q_out: (Vec<usize>, u32) = partial_search(&query.query, &fm_struct);
            query_arr[idx].hits = q_out.0;
            query_arr[idx].match_len = q_out.1;
        }
    } else {
        println!("Input query mode of: {} was invalid", query_mode);
        return 
    }

    // Output the query struct to output file
    // println!("{:?}", query_arr);
    println!("[*] Writing to file: {}", output);
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    for query in query_arr {
        if let Err(e) = writeln!(write_output, "{}\t{}\t{}\t{}", query.header, query.match_len, query.hits.len(),
                                 query.hits.iter().map( |&id| id.to_string() + "\t").collect::<String>()) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }
}

pub fn partial_search(query: &str, fm_struct: &FMindex) -> (Vec<usize>, u32) {
    let mut hits: Vec<usize> = Vec::new();
    let mut match_len: u32 = 1;
    let mut possible: (u32, u32) = (0, 0);
    let mut prev_possible: (u32, u32) = (0, 0);
    //A offset is just 1 since there's 1 '&'
    let c_offset: u32 = 1 + fm_struct.first_column[1];
    let g_offset: u32 = c_offset + fm_struct.first_column[2];
    let t_offset: u32 = g_offset + fm_struct.first_column[3];
    let rev_query: &str = &query.chars().rev().collect::<String>();
    let head_rest = split_first_char(rev_query).unwrap();

    //println!("Split into: {} and {}", head_rest.0, head_rest.1);
    // Grabbing index of rows in first column with head
    match head_rest.0 {
        'A' => possible = (1, c_offset),
        'C' => possible = (c_offset, g_offset),
        'G' => possible = (g_offset, t_offset),
        'T' => possible = (t_offset, fm_struct.bwt.len() as u32),
        _ => ()
    }

    // Iterate over current possibilities

    // Weird bug where my range is including one below hit
    // check for match of new char
    for c in head_rest.1.chars() {
        match_len += 1;
        prev_possible = possible;
        //println!("Now checking char: {}, possible size of: {:?}", c, possible.1 - possible.0);
        match c {
            'A' => possible = (occ(possible.0 as usize, 'A', &fm_struct.tally), occ((possible.1 - 1) as usize, 'A', &fm_struct.tally) + 1),
            'C' => possible = (c_offset + occ(possible.0 as usize, 'C', &fm_struct.tally) - 1, c_offset + occ((possible.1 - 1) as usize, 'C', &fm_struct.tally)),
            'G' => possible = (g_offset + occ(possible.0 as usize, 'G', &fm_struct.tally) - 1, g_offset + occ((possible.1 - 1) as usize, 'G', &fm_struct.tally)),
            'T' => possible = (t_offset + occ(possible.0 as usize, 'T', &fm_struct.tally) - 1, t_offset + occ((possible.1 - 1) as usize, 'T', &fm_struct.tally)),
            _ => ()
        }

        if possible.1 as i32 - possible.0 as i32 == 1 {
            match_len -= 1;
            for i in prev_possible.0 + 1..prev_possible.1 {
                hits.push(fm_struct.sa[i as usize]);
            }
            break;
        }

    }

    // Trivially returning match_len since in complete mode, its either query len or 0
    if possible.1 as i32 - possible.0 as i32 > 1 {
        match_len = query.len() as u32;
        for i in possible.0 + 1..possible.1 {
            hits.push(fm_struct.sa[i as usize]);
        }
    }

    return (hits, match_len);
}

pub fn complete_search(query: &str, fm_struct: &FMindex) -> (Vec<usize>, u32) {
    let mut hits: Vec<usize> = Vec::new();
    let mut match_len: u32 = 0;
    let mut possible: (u32, u32) = (0, 0);
    //A offset is just 1 since there's 1 '&'
    let c_offset: u32 = 1 + fm_struct.first_column[1];
    let g_offset: u32 = c_offset + fm_struct.first_column[2];
    let t_offset: u32 = g_offset + fm_struct.first_column[3];
    let rev_query: &str = &query.chars().rev().collect::<String>();
    let head_rest = split_first_char(rev_query).unwrap();

    //println!("Split into: {} and {}", head_rest.0, head_rest.1);
    // Grabbing index of rows in first column with head
    match head_rest.0 {
        'A' => possible = (1, c_offset),
        'C' => possible = (c_offset, g_offset),
        'G' => possible = (g_offset, t_offset),
        'T' => possible = (t_offset, fm_struct.bwt.len() as u32),
        _ => ()
    }

    // Iterate over current possibilities

    // Weird bug where my range is including one below hit
    // check for match of new char
    for c in head_rest.1.chars() {
        //println!("Now checking char: {}, possible size of: {:?}", c, possible.1 - possible.0);
        match c {
            'A' => possible = (occ(possible.0 as usize, 'A', &fm_struct.tally), occ((possible.1 - 1) as usize, 'A', &fm_struct.tally) + 1),
            'C' => possible = (c_offset + occ(possible.0 as usize, 'C', &fm_struct.tally) - 1, c_offset + occ((possible.1 - 1) as usize, 'C', &fm_struct.tally)),
            'G' => possible = (g_offset + occ(possible.0 as usize, 'G', &fm_struct.tally) - 1, g_offset + occ((possible.1 - 1) as usize, 'G', &fm_struct.tally)),
            'T' => possible = (t_offset + occ(possible.0 as usize, 'T', &fm_struct.tally) - 1, t_offset + occ((possible.1 - 1) as usize, 'T', &fm_struct.tally)),
            _ => ()
        }
    }

    // Trivially returning match_len since in complete mode, its either query len or 0
    if possible.1 as i32 - possible.0 as i32 > 1 {
        match_len = query.len() as u32;
        for i in possible.0 + 1..possible.1 {
            hits.push(fm_struct.sa[i as usize]);
        }
    }

    return (hits, match_len);
}

pub fn occ(row: usize, ch: char, tally: &Vec<Vec<u32>>) -> u32 {
    match ch {
        'A' => return tally[row][1],
        'C' => return tally[row][2],
        'G' => return tally[row][3],
        'T' => return tally[row][4],
        _ => 0
    }
}

// Splitting string into head and rest
pub fn split_first_char(s: &str) -> Option<(char, &str)> {
    let mut chars = s.chars();
    chars.next().map(|c| (c, chars.as_str()))
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
pub fn read_queries(path: &str) -> Vec<FMQ> {
    let query_file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(query_file);
    let mut query_arr: Vec<FMQ> = Vec::new();
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
                            FMQ {
                                header: header.to_string(),
                                query: query.to_string(),
                                hits: Vec::new(),
                                match_len: 0,
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
        FMQ {
            header: header.to_string(),
            query: query.to_string(),
            hits: Vec::new(),
            match_len: 0,
        }
    });
    return query_arr;
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
