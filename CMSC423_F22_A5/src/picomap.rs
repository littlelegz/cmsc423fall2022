use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;
use serde::{Serialize, Deserialize};
use bincode;
use std::cmp::{Ordering, min, max};
use priority_queue::PriorityQueue;
use std::collections::HashMap;

// Structs
#[derive(Debug, Clone)]
pub struct MapQuery {
    problem_name: String,
    query: String,
    score: i32,
    ref_start: Vec<usize>,
    cigar: Vec<String>
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
        // Finding seeds
        let exact_match: (Vec<usize>, u32) = partial_search(&query.query, &reference);

        // Get reference slice

        // Calling aligner
        if exact_match.1 == 100 {
            let count = exact_match.0.len();
            problems[idx].ref_start = exact_match.0;
            problems[idx].score = 0;
            for i in 0..count {
                problems[idx].cigar.push("100=".to_string());
            }
        } else {// No exact match found
            let mut fmseed: Vec<usize> = find_seeds(&query.query, &reference);
            //println!("Found fmseeds: {:?}", fmseed);
            let mut best_seed: usize = 0;
            let mut best_score: i32 = -1000000;
            let mut result_cigar: Vec<String> = Vec::new();
            let mut more_refs: Vec<usize> = Vec::new();
            for seed in fmseed {
                let lower: usize;
                if seed > 15 { // accounting for the 15 bases before the seed
                    lower = seed - 15;
                } else {
                    lower = 0;
                }
                let reference_slice: &str = &reference.genome[lower..seed+query.query.len()+15];
                let result: (i32, usize, usize, String) = reference_align(query, reference_slice, -mismatch_pen, -gap_pen);
                if result.0 > best_score {
                    best_score = result.0;
                    best_seed = seed - 15 + result.1;
                    result_cigar = vec![result.3];
                    more_refs = vec![best_seed];
                } else if result.0 == best_score {
                    more_refs.push(seed - 15 + result.1);
                    result_cigar.push(result.3);
                }
            }
            
            /*let curr_start = best_seed;
            let new_ref: &str = &reference.genome[curr_start..curr_start+query.query.len()];
            let more_refs: (Vec<usize>, u32) = partial_search(new_ref, &reference);

            if more_refs.1 == query.query.len() as u32{
                problems[idx].ref_start = more_refs.0;
                problems[idx].score = best_score;
                problems[idx].cigar = result_cigar;
            } else {
                problems[idx].ref_start = vec![best_seed];
                problems[idx].score = best_score;
                problems[idx].cigar = result_cigar;
            }*/
            more_refs.dedup();
            problems[idx].ref_start = more_refs;
            problems[idx].score = best_score;
            problems[idx].cigar = result_cigar;
        }
    }

    //let mut cigar: Vec<(&str, &str)> = Vec::new(); 

    //println!("Problem: {:?}", problem);

    // Creating output file
    println!("[*] Writing to file: {}", output);
    if Path::new(&output).exists() {
        fs::remove_file(&output).unwrap();
    }
    let mut write_output = OpenOptions::new().append(true).create(true).open(output).expect("Unable to open file"); 
    for query in problems {
        if let Err(e) = writeln!(write_output, "{}\t{}", query.problem_name, query.ref_start.len()) {
            eprintln!("Couldn't write to file: {}", e);
        }
        for (index, i) in query.ref_start.iter().enumerate() {
            if let Err(e) = writeln!(write_output, "{}\t{}\t{}", i, query.score, query.cigar[index as usize]) {
                eprintln!("Couldn't write to file: {}", e);
            }
        }
    }

}

// Given the query and fm index, return a list of possible seeds, sorted by most likely to least likely

pub fn find_seeds(query: &str, fm_struct: &FMindex) -> Vec<usize> {
    // Divide string into fouths
    let mut seeds: Vec<(usize, u32)> = Vec::new();
    let mut occ_seeds: HashMap<usize, u32> = HashMap::new();
    let seed_len: usize = query.len() / 4;
    //println!("Query len: {}, Seed length: {}", query.len(), seed_len);
    let mut seed: &str = &query[0..seed_len];
    // 1/4 of the query
    let mut fmseed: (Vec<usize>, u32) = partial_search(seed, fm_struct);
    // Before adjusting, must remove seeds which result in subtract overflow
    fmseed.0.retain(|i| {
        let delete = {
            if *i as i32 - (seed_len as i32 - fmseed.1 as i32) < 0 {
                true
            } else {
                false
            }
        };
        !delete
    });
    let mut adjustedseeds: Vec<usize> = fmseed.0.iter().map(|&x| x - (seed_len - fmseed.1 as usize)).collect();
    let mut adjust_weight_seeds: Vec<(usize, u32)> = adjustedseeds.iter().map(|&x| (x, fmseed.1)).collect();
    seeds.append(&mut adjust_weight_seeds);
    //println!("Found seeds 1/4: {:?}", seeds);
    // 2/4 of the query
    seed = &query[seed_len..seed_len * 2];
    fmseed = partial_search(seed, fm_struct);
    fmseed.0.retain(|i| {
        let delete = {
            if (*i as i32 - (seed_len as i32 - fmseed.1 as i32) - seed_len as i32) < 0 {
                true
            } else {
                false
            }
        };
        !delete
    });
    adjustedseeds = fmseed.0.iter().map(|&x| x - (seed_len - fmseed.1 as usize) - seed_len).collect();
    adjust_weight_seeds = adjustedseeds.iter().map(|&x| (x, fmseed.1)).collect();
    seeds.append(&mut adjust_weight_seeds);
    //println!("Found seeds 2/4: {:?}", seeds);
    // 3/4 of the query
    seed = &query[seed_len * 2..seed_len * 3];
    fmseed = partial_search(seed, fm_struct);
    fmseed.0.retain(|i| {
        let delete = {
            if (*i as i32 - (seed_len as i32 - fmseed.1 as i32) - (seed_len*2) as i32) < 0 {
                true
            } else {
                false
            }
        };
        !delete
    });
    adjustedseeds = fmseed.0.iter().map(|&x| x - (seed_len - fmseed.1 as usize) - seed_len*2).collect();
    adjust_weight_seeds = adjustedseeds.iter().map(|&x| (x, fmseed.1)).collect();
    seeds.append(&mut adjust_weight_seeds);
    //println!("Found seeds 3/4: {:?}", seeds);

    // 4/4 of the query
    seed = &query[seed_len * 3..];
    fmseed = partial_search(seed, fm_struct);
    fmseed.0.retain(|i| {
        let delete = {
            if (*i as i32 - (seed_len as i32 - fmseed.1 as i32) - (seed_len*3) as i32) < 0 {
                true
            } else {
                false
            }
        };
        !delete
    });
    //println!("Found unadjusted seeds 4/4: {:?}", fmseed.0);
    adjustedseeds = fmseed.0.iter().map(|&x| x - (seed_len - fmseed.1 as usize) - seed_len*3).collect();
    adjust_weight_seeds = adjustedseeds.iter().map(|&x| (x, fmseed.1)).collect();
    seeds.append(&mut adjust_weight_seeds);
    //println!("Found seeds 4/4: {:?}", seeds);

    // insert all seeds into an occurance hashmap
    for seed in seeds {
        if occ_seeds.contains_key(&seed.0) {
            *occ_seeds.get_mut(&seed.0).unwrap() += seed.1;
        } else {
            occ_seeds.insert(seed.0, seed.1);
        }
    }
    // Use priority queue to sort by most occuring
    let mut pq: PriorityQueue<usize, u32> = PriorityQueue::new();
    for key in occ_seeds.keys() {
        pq.push(*key, occ_seeds[key]);
    }
    let top_priority: u32 = *pq.peek().unwrap().1;
    let mut output: Vec<usize> = Vec::new();
    while pq.peek() != None && *pq.peek().unwrap().1 == top_priority {
        output.push(pq.pop().unwrap().0);
    }

    return output
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
    let mut cigar: String = String::new();

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
    let mut query: String = String::new();
    let mut header: String = String::new();
    
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
                ref_start: Vec::new(),
                cigar: Vec::new()
            });
            //Reset header and query
            header = String::new();
            query = String::new();
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
