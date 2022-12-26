use std::env;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use json::JsonValue;

fn main() {
    let args: Vec<String> = env::args().collect();
    let path = &args[1];
    let file = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(file);
    let mut firstrun = true;
    let mut tot_len = 0;
    let mut min_len = u32::MAX;
    let mut max_len = 0;
    let mut num_records = 0;
    let mut agg_length = 0;
    let mut stats = json::parse(r#"
        {
            "min_len": 0,
            "max_len": 0,
            "mean_len": 0,
            "tot_len": 0,
            "num_records": 0
        }
    "#).unwrap();

    for line in reader.lines() {
        if let Ok(curr) = line {
            match curr.chars().next() {
                Some(firstchar) => {
                    if firstchar == '>' {
                        if !firstrun {
                            if agg_length > max_len {
                                max_len = agg_length;
                            }
                            if agg_length < min_len {
                                min_len = agg_length;
                            }
                        }
                        num_records += 1;
                        agg_length = 0;
                    } else {
                        firstrun = false;
                        tot_len += curr.chars().count();
                        agg_length = agg_length + curr.len() as u32;
                    }
                }
                None => {
                    break;
                }
            }
            
        }
    }
    // Checking stats on final read

    if agg_length > max_len {
        max_len = agg_length;
    }
    if agg_length < min_len {
        min_len = agg_length;
    }

    stats["tot_len"] = JsonValue::from(tot_len);
    stats["max_len"] = JsonValue::from(max_len);
    stats["min_len"] = JsonValue::from(min_len);
    stats["num_records"] = JsonValue::from(num_records);
    if num_records != 0 {
        stats["mean_len"] = JsonValue::from(tot_len/num_records);
    }
    println!("{:#}", stats);
}
