use std::env;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::process::exit;

#[derive(Debug, Clone)]
pub struct QP {
    problem_name: String,
    query_x: String,
    query_y: String,
    score: i32,
    y_start: usize,
    y_end: usize,
    cigar: String
}


fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() != 6 {
        println!("Invalid parameters");
        exit(1);
    }

    let mut problems: Vec<QP> = read_queries(&args[1]);
    let method: &str = &args[2];
    let mismatch_pen: i32 = args[3].parse().unwrap();
    let gap_pen: i32 = args[4].parse().unwrap();
    let output: &str = &args[5];

    if method == "global" {
        for (idx, query) in problems.clone().iter().enumerate() {
            let result: (i32, usize, usize, String) = global_align(query, -mismatch_pen, -gap_pen);
            problems[idx].score = result.0;
            problems[idx].y_start = result.1;
            problems[idx].y_end = result.2;
            problems[idx].cigar = result.3;
        }
    } else if method == "fitting" {
        for (idx, query) in problems.clone().iter().enumerate() {
            let result: (i32, usize, usize, String) = fitting_align(query, -mismatch_pen, -gap_pen);
            problems[idx].score = result.0;
            problems[idx].y_start = result.1;
            problems[idx].y_end = result.2;
            problems[idx].cigar = result.3;
        }
    } else {
        println!("Invalid method");
        exit(1);
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
        if let Err(e) = writeln!(write_output, "{}\n{}\n{}\n{}\t{}\t{}\t{}" , query.problem_name, query.query_x, query.query_y, query.score, query.y_start, query.y_end, query.cigar) {
            eprintln!("Couldn't write to file: {}", e);
        }
    }

}

pub fn global_align(query: &QP, mismatch_pen: i32, gap_pen: i32) -> (i32, usize, usize, String) {
    // convert queries to a vector of chars
    let x: Vec<char> = query.query_x.chars().collect();
    let y: Vec<char> = query.query_y.chars().collect();

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

    // initialize first row and column
    for i in 0..m + 1{
        matrix[i][0] = (i as i32 * gap_pen, 0);
    }
    for j in 0..n + 1{
        matrix[0][j] = (j as i32 * gap_pen, 0);
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

    //println!("{:?}", matrix[m][n]);
    
    // constructing the cigar string
    let mut i: usize = m;
    let mut j: usize = n;
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
    } else if j > 0 {
        for _ in 0..j {
            cigar = "D".to_string() + &cigar;
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

    return (matrix[m][n].0, 0, n, cigar_str);
}

pub fn fitting_align(query: &QP, mismatch_pen: i32, gap_pen: i32) -> (i32, usize, usize, String) {
    // convert queries to a vector of chars
    let x: Vec<char> = query.query_x.chars().collect();
    let y: Vec<char> = query.query_y.chars().collect();

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

pub fn read_queries(path: &str) -> Vec<QP> {
    let query_file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(query_file);
    let mut query_arr: Vec<QP> = Vec::new();
    let mut query: String = "".to_owned();
    let mut header: String = "".to_owned();
    
    for line in reader.lines() {
        if header == "" {
            header = line.unwrap(); //No header? Read header
        } else if query == "" {
            query = line.unwrap(); //No query? Read query
        } else {
            query_arr.push(QP { //Have header and query, push to array
                problem_name: header,
                query_x: query,
                query_y: line.unwrap(),
                score: 0,
                y_start: 0,
                y_end: 0,
                cigar: "".to_string()
            });
            //Reset header and query
            header = "".to_owned();
            query = "".to_owned();
        }
    }

    return query_arr;
}