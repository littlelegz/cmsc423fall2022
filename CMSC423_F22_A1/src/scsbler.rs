use std::env;
use std::process::exit;
use std::fs::{File, OpenOptions};
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use std::collections::HashMap;
use std::iter;
use std::path::Path;
use priority_queue::PriorityQueue;
use std::cmp;

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() != 4 {
        println!("Invalid parameters");
        exit(1);
    }
    let reads = &args[1];
    let min_olap: u32 = args[2].parse().unwrap();
    let output_stem = &args[3];
    // Creating output files
    let output_fa = format!("{}.fa", output_stem);
    let output_edges = format!("{}.edges", output_stem);
    if Path::new(&output_fa).exists() {
        fs::remove_file(&output_fa).unwrap();
    }
    if Path::new(&output_edges).exists() {
        fs::remove_file(&output_edges).unwrap();
    }
    let mut write_fa = OpenOptions::new().append(true).create(true).open(output_fa).expect("Unable to open file");   
    let mut write_edges = OpenOptions::new().append(true).create(true).open(output_edges).expect("Unable to open file");   

    // Read reads into reads vec
    let genome_file = File::open(reads).expect("Unable to open file");
    let mut reads: Vec<String> = Vec::new();
    let reader = BufReader::new(genome_file);
    let mut insert_line = "".to_string();
    for line in reader.lines() {
        if let Ok(curr) = line {
            match curr.chars().next() {
                Some(firstchar) => {
                    if firstchar == '>' {
                        if !insert_line.is_empty() {
                            reads.push(insert_line);
                        }
                        insert_line = "".to_string();
                        continue;
                    } else {
                        insert_line = format!("{}{}", insert_line, curr);
                    }
                }
                None => {
                    break;
                }
            }
            
        }
    }
    // Add final line
    reads.push(insert_line);
    println!("[*] Read {} lines from provided filepath into reads", reads.len());
    println!("[*] Computing overlap graph (this takes some time)...");

    // Graph datastructure
    // let mut overlap_graph: HashMap<i32, Vec<Vec<i32>>> = HashMap::new();
    // Union find datastructure
    let mut uni_find: HashMap<u32, Option<u32>> = HashMap::new();
    // Suffix/Prefix hash
    let mut fix_hash: HashMap<String, Vec<u32>> = HashMap::new();
    // Priority queue to store overlaps by rank
    let mut ranked_overlaps: PriorityQueue<Vec<u32>, Vec<i32>> = PriorityQueue::new();
    // Calculate overlaps
    // Storing prefixes and suffixes as key with value being read id
    for i in 0..reads.len() {
        // Insert read ids into uni_find
        uni_find.insert(i as u32, None);
        let curr = &reads[i];
        let prefix = &curr[..min_olap as usize];
        let suffix = &curr[(curr.chars().count() - min_olap as usize)..];
        if fix_hash.get_mut(prefix) == None {
            fix_hash.insert(prefix.to_string(), vec![i as u32]);
        } else {
            fix_hash.get_mut(prefix).unwrap().push(i as u32);
        }
        if fix_hash.get_mut(suffix) == None {
            fix_hash.insert(suffix.to_string(), vec![i as u32]);
        } else {
            fix_hash.get_mut(suffix).unwrap().push(i as u32);
        }
    }
    // Using l-mers to find overlaps
    for i in 0..reads.len() {
        let mut l_start = 0;
        let curr = &reads[i];
        // Lookup l-mer
        while l_start + min_olap <= curr.chars().count() as u32 {
            let l_mer = &curr[l_start as usize..(l_start + min_olap) as usize];
            if fix_hash.get_mut(l_mer) != None {
                let v: &mut Vec<u32> = fix_hash.get_mut(l_mer).unwrap();
                for v_i in v {
                    if *v_i == i as u32 {
                        continue;
                    }
                    let overlap = find_overlap(curr, &reads[*v_i as usize], min_olap);
                    // insert into pq
                    if overlap >= min_olap {
                        ranked_overlaps.push(vec![i as u32, *v_i], vec![overlap as i32, -(i as i32), -(*v_i as i32)]); 
                    }
                }
            }
            l_start += 1;
        }
    }
    
    println!("Graph computed and written to file");
    // Peek into the priority queue
    println!("First item in priority queue is: {:?}", ranked_overlaps.peek());
    // Work through all pq entries
    let pq = ranked_overlaps.clone();
    println!("[*] Starting priority queue: ");
    for x in pq.into_sorted_iter() {
        if let Err(e) = writeln!(write_edges, "{}\t{}\t{}", x.0[0], x.0[1], x.1[0]) {
            eprintln!("Couldn't write to file: {}", e);
        }
        //println!("{:?}", x)
    } 

    while ranked_overlaps.peek() != None {
        let curr = ranked_overlaps.pop().unwrap();
        let mut from = curr.0[0];
        let mut to = curr.0[1];
        println!("Popping element: {:?}", curr);
        // Combine strings and save to reads
        // Case: Two fresh strings
        if uni_find.get(&from) == None && uni_find.get(&to) == None {
            let suffix = &reads[from as usize];
            let prefix = &reads[to as usize][(curr.1[0] as usize)..];
            let extended = format!("{}{}", suffix, prefix);
            let new_pos = cmp::min(from, to);
            let old_pos = cmp::max(from, to);
            uni_find.insert(old_pos, Some(new_pos));
            // Update combined into reads with smallest rank
            reads[new_pos as usize] = extended;
            // Delete j read since it will now be referenced by i read
            // Attempt to remove the back edge if it exists
            ranked_overlaps.remove(&vec![curr.0[1], curr.0[0]]);
            reads[old_pos as usize] = "".to_string();
        } else { // Case either string has already been merged
            while uni_find.get(&from) != None {
                if *(uni_find.get(&from).unwrap()) == None {
                    break;
                }
                from = uni_find.get(&from).unwrap().unwrap();
            }
            while uni_find.get(&to) != None {
                if *(uni_find.get(&to).unwrap()) == None {
                    break;
                }
                to = uni_find.get(&to).unwrap().unwrap();
            }
            if from == to {
                continue;
            }
            let check_overlap = find_overlap(&reads[from as usize], &reads[to as usize], min_olap);
            if curr.1[0] <= check_overlap as i32 {
                let suffix = &reads[from as usize];
                let prefix = &reads[to as usize][(curr.1[0] as usize)..];
                let extended = format!("{}{}", suffix, prefix);
                let new_pos = cmp::min(from, to);
                let old_pos = cmp::max(from, to);
                uni_find.insert(old_pos, Some(new_pos));
                // Update combined into reads with smallest rank
                reads[new_pos as usize] = extended;
                // Delete j read since it will now be referenced by i read
                // Attempt to remove the back edge if it exists
                ranked_overlaps.remove(&vec![curr.0[1], curr.0[0]]);
                reads[old_pos as usize] = "".to_string();
            } else {
                if check_overlap >= min_olap {
                    ranked_overlaps.push(vec![from as u32, to as u32], vec![check_overlap as i32, -(from as i32), -(to as i32)]);
                }
            }
        }

        // let pq = ranked_overlaps.clone();
        /* println!("[*] Current priority queue: ");
        for x in pq.into_sorted_iter() {
            println!("{:?}", x)
        }  */
    }

    let mut island_count = 0;
    for x in 0..reads.len() {
        if !&reads[x].is_empty() {
            if let Err(e) = writeln!(write_fa, ">{}:{}", island_count, reads[x].chars().count()) {
                eprintln!("Couldn't write to file: {}", e);
            }
            island_count += 1;
            if let Err(e) = writeln!(write_fa, "{}", &reads[x]) {
                eprintln!("Couldn't write to file: {}", e);
            }
            //println!("Generated SCS read at {} with len {}", x, &reads[x]);
        }
    }
    println!("Generated islands: {}", island_count);

}

pub fn find_overlap(curr: &str, comp: &str, min_olap: u32) -> u32 {
    let curr_suffixes = suffixes(&curr).collect::<Vec<_>>();
    let comp_prefixes = prefixes(&comp).collect::<Vec<_>>();
    let max_overlap: u32;
    if curr_suffixes.len() > comp_prefixes.len() {
        max_overlap = comp_prefixes.len() as u32;
    } else {
        max_overlap = curr_suffixes.len() as u32;
    }
            
    for x in (min_olap..max_overlap).rev() {
        if curr_suffixes[x as usize].eq(comp_prefixes[x as usize]) {
            return x;
        }
    }
    return 0;
}

pub fn prefixes(s: &str) -> impl Iterator<Item = &str> + DoubleEndedIterator {
    s.char_indices()
        .map(move |(pos, _)| &s[..pos])
        .chain(iter::once(s))
}

pub fn suffixes(s: &str) -> impl Iterator<Item = &str> + DoubleEndedIterator {
    s.char_indices()
        .map(move |(pos, _)| &s[pos..])
        .chain(iter::once(""))
        .rev()
}