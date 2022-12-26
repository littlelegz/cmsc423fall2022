use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;
use std::cmp::{Ordering, min, max};

// Structs
#[derive(Debug, Clone)]
pub struct MapQuery {
    problem_name: String,
    query: String,
    score: i32,
    ref_start: usize,
    cigar: String
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FMindex {
    genome: String,
    bwt: Vec<String>,
    sa: Vec<usize>,
    first_column: Vec<u32>,
    tally: Vec<Vec<u32>>,
}


fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() != 6 {
        println!("Invalid parameters");
        exit(1);
    }

    // Read input args
    let reference: FMindex = read_fm_struct(&args[1]);
    let mut problems: Vec<MapQuery> = read_queries(&args[2]);
    let mismatch_pen: i32 = args[3].parse().unwrap();
    let gap_pen: i32 = args[4].parse().unwrap();
    let output: &str = &args[5];

    // Map queries
    for (idx, query) in problems.clone().iter().enumerate() {
        // Find seed (best partial mapping)
        //println!("[*] Now calculating full alignment for query {}", query.problem_name);
        let fmseed: (Vec<usize>, u32) = partial_search(&query.query, &reference);
        let saseed: (Vec<usize>, u32, u32);
        //println!("Found fmseeds: {:?}", fmseed);
        // Convert seed into string of reference
        // Takes second found if first seed index out of range
        let mut reference_slice: &str = "";

        /*if (fmseed.0[0] as i32 - (query.query.len() - fmseed.1 as usize) as i32) < 0 {
            fmseed.0.remove(0);
        }*/

        if fmseed.1 as usize == query.query.len() { // Found exact match, no need for align
            problems[idx].ref_start = fmseed.0[0];
            problems[idx].score = 0;
            problems[idx].cigar = format!("{}{}", query.query.len(), "="); 
        } else if fmseed.0.len() > 0 && fmseed.0[0] as i32 - (query.query.len() as i32 - fmseed.1 as usize as i32) > 0 { // Found partial match, need to align
            // Setting reference to slice of reference genome
            reference_slice = &reference.genome[fmseed.0[0] - (query.query.len() - fmseed.1 as usize)..
                                            fmseed.0[0] + fmseed.1 as usize + 15];
            problems[idx].ref_start = fmseed.0[0] - (query.query.len() - fmseed.1 as usize);
        
        // Calling aligner
            let mut result: (i32, usize, usize, String) = reference_align(query, reference_slice, -mismatch_pen, -gap_pen);
            problems[idx].score = result.0;
            problems[idx].cigar = result.3;
        } else { // No proper seeds from fm search, using saseeds
            saseed = binary_search_wrapper(&query.query, reference.sa.clone(), &reference.genome);
            reference_slice = &reference.genome[saseed.0[0]..saseed.0[0] + query.query.len()];
            let result: (i32, usize, usize, String) = reference_align(query, reference_slice, -mismatch_pen, -gap_pen);
            problems[idx].ref_start = saseed.0[0];
            problems[idx].score = result.0;
            problems[idx].cigar = result.3;
        }
        
    }

    // Creating output file
    println!("[*] Writing to file: {}", output);
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    for query in problems {
        if let Err(e) = writeln!(write_output, "{}\t1\n{}\t{}\t{}", query.problem_name, query.ref_start, query.score, query.cigar) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }

}
// Binary search into sa
pub fn binary_search_wrapper(query: &str, sa: Vec<usize>, genome: &str) -> (Vec<usize>, u32, u32) {
    // Perform lower search
    let lower = binary_search(query, sa.clone(), genome);
    let _lower_string: &str = &genome[sa[lower.0 as usize]..];
    //println!("[*] Found lower at: {}, Upper at: {}, Sa lower: {:?}, Sa Upper: {:?}", lower.0, upper.0, sa[lower.0 as usize], sa[upper.0 as usize]);
    
    return (vec![sa[lower.0 as usize]], lower.1, 0);
}

pub fn binary_search(target: &str, sa: Vec<usize>, genome: &str) -> (u32, u32) {
    let mut low: i32 = 0;
    let mut high: i32 = sa.len() as i32;
    let mut comparisons: u32 = 0;
    let mut lcp: u32 = 0;
    let mut best: u32 = 0;
 
    loop {
        let middle: i32 = (high + low) / 2;
        let current = &genome[sa[middle as usize]..];
        let result = compare(target, &current);
        if result.1 - 1 > lcp {
            lcp = result.1 - 1;
            best = middle as u32;
        }
        comparisons += result.1;

        if result.0 == Ordering::Less {
            if middle == low + 1 {
                return (best, comparisons);
            }
            high = middle
        }
        if result.0 == Ordering::Greater {
            if middle == high - 1 {
                return (best, comparisons);
            }
            low = middle
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

pub fn reference_align(query: &MapQuery, ref_start: &str, mismatch_pen: i32, gap_pen: i32) -> (i32, usize, usize, String) {
    // convert queries to a vector of chars
    let x: Vec<char> = query.query.chars().collect();
    let y: Vec<char> = ref_start.chars().collect();

    // store length of x and y
    let m: usize = x.len();
    let n: usize = y.len();

    // initialize alignment matrix
    let mut matrix: Vec<Vec<(i32, i32)>> = Vec::new();
    for _ in 0..m + 1 {
        let mut row: Vec<(i32, i32)> = Vec::new();
        for _ in 0..n + 1 {
            row.push((0, 0));
        }
        matrix.push(row);
    }

    // ignore gaps before X
    for i in 0..m + 1{
        matrix[i][0] = (i as i32 * gap_pen, 0);
    }
    for j in 0..n + 1{
        matrix[0][j] = (0, 0);
    }

    //fill matrix 1 is diag, 2 is up, 3 is left
    for i in 1..m + 1 {
        for j in 1..n + 1 {
            let mut max: i32 = 0;
            let mut max_idx: i32 = 0;
            let mut score: i32 = 0;
            if x[i-1] == y[j-1] {
                score = 0;
            } else {
                score = mismatch_pen;
            }
            let diag: i32 = matrix[i - 1][j - 1].0 + score;
            let left: i32 = matrix[i - 1][j].0 + gap_pen;
            let up: i32 = matrix[i][j - 1].0 + gap_pen;
            if diag >= left && diag >= up {
                max = diag;
                max_idx = 1;
            } else if left >= diag && left >= up {
                max = left;
                max_idx = 2;
            } else {
                max = up;
                max_idx = 3;
            }
            matrix[i][j] = (max, max_idx);
        }
    }

    // find traceback start point
    let mut max: i32 = -999999999;
    let mut max_i: usize = 0;
    let mut max_j: usize = 0;
    for j in 1..n + 1 {
        if matrix[m][j].0 > max {
            max = matrix[m][j].0;
            max_i = m;
            max_j = j;
        }
    }

    //println!("Found max of {} at ({}, {})", max, max_i, max_j);
    let mut i: usize = max_i;
    let mut j: usize = max_j;
    let mut cigar: String = "".to_string();

    while i > 0 && j > 0 {
        if matrix[i][j].1 == 1 {
            if x[i-1] == y[j-1] { // diag
                cigar = "=".to_string() + &cigar;
            } else {
                cigar = "X".to_string() + &cigar;
            }
            i -= 1;
            j -= 1;
        } else if matrix[i][j].1 == 2 { // left
            cigar = "I".to_string() + &cigar;
            i -= 1;
        } else { // up
            cigar = "D".to_string() + &cigar;
            j -= 1;
        }
    }

    // Add remaining gaps
    if i > 0 {
        for _ in 0..i {
            cigar = "I".to_string() + &cigar;
        }
    }
    // println!("{:?}", cigar);
    // convert cigar string into a vector of tuples
    let mut cigar_vec: Vec<(&str, String)> = Vec::new();
    let mut count: usize = 0;
    let mut prev: &str = &cigar[0..1];
    for i in 0..cigar.len() {
        if &cigar[i..i+1] == prev {
            count += 1;
        } else {
            cigar_vec.push((prev, count.to_string()));
            prev = &cigar[i..i+1];
            count = 1;
        }
    }
    // push last element
    cigar_vec.push((prev, count.to_string()));
    //println!("{:?}", cigar_vec);

    // convert cigar vector into a string
    let mut cigar_str: String = "".to_string();
    for (c, n) in cigar_vec {
        cigar_str = cigar_str + &n + c;
    }


    return (matrix[max_i][max_j].0, j, max_j, cigar_str);
}

pub fn read_queries(path: &str) -> Vec<MapQuery> {
    let query_file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(query_file);
    let mut query_arr: Vec<MapQuery> = Vec::new();
    let mut query: String = "".to_owned();
    let mut header: String = "".to_owned();
    
    for line in reader.lines() {
        if header == "" {
            header = line.unwrap(); //No header? Read header
            header.remove(0); 
        } else {
            query = line.unwrap(); //No query? Read query and push to array
            query_arr.push(MapQuery { //Have header and query, push to array
                problem_name: header,
                query: query,
                score: 0,
                ref_start: 0,
                cigar: "".to_string()
            });
            //Reset header and query
            header = "".to_owned();
            query = "".to_owned();
        }
    }

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

pub fn occ(row: usize, ch: char, tally: &Vec<Vec<u32>>) -> u32 {
    match ch {
        'A' => return tally[row][1],
        'C' => return tally[row][2],
        'G' => return tally[row][3],
        'T' => return tally[row][4],
        _ => 0
    }
}

pub fn split_first_char(s: &str) -> Option<(char, &str)> {
    let mut chars = s.chars();
    chars.next().map(|c| (c, chars.as_str()))
}

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
