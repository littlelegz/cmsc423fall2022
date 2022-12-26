use std::env;
use std::fs::{File, OpenOptions};
use std::fs;
use std::io::{prelude::*, BufReader, Write};
use json::JsonValue;
use std::process::exit;
use rand::Rng;
use std::path::Path;

fn main() {
    // Reading args
    let args: Vec<_> = env::args().collect();
    if args.len() != 6 {
        println!("Invalid parameters");
        exit(1);
    }
    // Reading inputs
    let read_len: u32 = args[1].parse().unwrap();
    let target_depth: u32 = args[2].parse().unwrap();
    let theta: f32 = args[3].parse().unwrap();
    let genome = &args[4];
    let output_stem = &args[5];
    // Creating output files
    let output_fa = format!("{}.fa", output_stem);
    let output_stats = format!("{}.stats", output_stem);
    if Path::new(&output_fa).exists() {
        fs::remove_file(&output_fa).unwrap();
    }
    if Path::new(&output_stats).exists() {
        fs::remove_file(&output_stats).unwrap();
    }
    let mut write_fa = OpenOptions::new().append(true).create(true).open(output_fa).expect("Unable to open file");   
    let mut write_stats = OpenOptions::new().append(true).create(true).open(output_stats).expect("Unable to open file");   
    // Reading genome into gen
    let genome_file = File::open(genome).expect("Unable to open file");
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
    println!("[*] Read genome into gen");
    let gen_len = gen.chars().count() as u32;
    println!("[*] Size of gen {}", gen_len);
    // Calculating m (required # of reads to achieve target_depth)
    let mut gen_array: Vec<i32> = vec![0; gen_len as usize];
    let mut m: f32 = ((target_depth * gen_len) / read_len) as f32;
    m = m.ceil();

    println!("[*] Desired m: {}", m);
    println!("[*] Params are:\n  > read_len: {}\n  > target_depth: {}\n  > theta: {}\n  > genome: {}\n  > output_stem: {}\n", read_len, target_depth, theta, genome, output_stem);
    
    // Generating reads and writing header + read to output_fa
    // Calculating bases_covered using ranges

    let mut ranges: Vec<Vec<u32>> = Vec::new();
    for x in 1..=(m as u32) {
        let num = rand::thread_rng().gen_range(0..(gen_len - read_len));
        let header = format!(">{}:{}:{}", x-1, num, read_len);
        let slice = &gen[(num as usize)..((num + read_len) as usize)];
        for i in num..num+read_len {
            gen_array[i as usize] += 1;
        }
        let range: Vec<u32> = vec![num, num+read_len];
        let mut to_add: bool = true;
        // Push ranges if theta overlap
        for i in 0..ranges.len() {
            if range[0] <= ranges[i][0] && range[1] >= (ranges[i][0] + (read_len as f32 *theta) as u32) && range[1] <= ranges[i][1] {
                to_add = false;
                ranges[i][0] = range[0];
            } else if range[1] >= ranges[i][1] && range[0] >= ranges[i][0] && range[0] <= (ranges[i][1] - (read_len as f32 *theta) as u32) {
                to_add = false;
                ranges[i][1] = range[1];
            } else if range[0] >= ranges[i][0] && range[1] <= ranges[i][1] {
                to_add = false;
            }
        }
        if to_add {
            ranges.push(range);
        }
        // Remove unnecessary ranges
        let mut curr = 0;
        let runs = ranges.len();
        for _a in 0..runs {
            let mut to_remove = false;
            for i in 0..ranges.len() {
                if i != curr {
                    if ranges[curr][0] <= ranges[i][0] && ranges[curr][1] >= (ranges[i][0] + (read_len as f32 *theta) as u32) && ranges[curr][1] <= ranges[i][1]  {
                        to_remove = true;
                        ranges[i][0] = ranges[curr][0];
                    } else if ranges[curr][1] >= ranges[i][1] && ranges[curr][0] >= ranges[i][0] && ranges[curr][0] <= (ranges[i][1] - (read_len as f32 *theta) as u32) {
                        to_remove = true;
                        ranges[i][1] = ranges[curr][1];
                    } else if ranges[curr][0] >= ranges[i][0] && ranges[curr][1] <= ranges[i][1] {
                        to_remove = true;
                    }
                }
            }
            if to_remove {
                ranges.remove(curr);
            } else {
                curr += 1;
            }
            
        }
        //println!("Rand: {}", num);
        if let Err(e) = writeln!(write_fa, "{}", header) {
            eprintln!("Couldn't write to file: {}", e);
        }
        if let Err(e) = writeln!(write_fa, "{}", slice) {
            eprintln!("Couldn't write to file: {}", e);
        }

    }

    
    println!("[*] Generated ranges: {:?}", ranges);

    // Generating stats json output
    let mut stats = json::parse(r#"
        {
            "num_reads": 0,
            "bases_covered": 0,
            "avg_depth": 0,
            "var_depth": 0,
            "num_islands": 0
        }
    "#).unwrap();
    stats["num_reads"] = JsonValue::from(m as u32);
    let avg_depth: f32 = m as f32 * read_len as f32 / gen_len as f32;
    stats["avg_depth"] = JsonValue::from(avg_depth);
    println!("[*] Calculated avg_depth: {:.5}", stats["avg_depth"]);

    // Calculate bases_covered
    let mut covered_count = 0;
    for (_pos, e) in gen_array.iter().enumerate() {
        if *e != 0 {
            covered_count += 1;
        }
    }
    println!("[*] Calculated coverage: {}", covered_count as f32 /gen_len as f32);
    // Calculate var_depth
    let mut var: f32 = 0.0;
    for (_pos, e) in gen_array.iter().enumerate() {
        var += (*e as f32 - avg_depth as f32).powf(2.0);
    }
    println!("[*] Calculated var: {}", var);
    let var_depth: f32 = var as f32 / (gen_len as f32 - 1.0) as f32;
    println!("[*] Calculated Var_Depth: {:.5}", var_depth);

    stats["bases_covered"] = JsonValue::from(covered_count);
    stats["var_depth"] = JsonValue::from(var_depth);
    stats["num_islands"] = JsonValue::from(ranges.len());

    if let Err(e) = writeln!(write_stats, "{:#}", stats) {
        eprintln!("Couldn't write to file: {}", e);
    }
}
