use std::thread;
use crossbeam_channel::{bounded, Sender, Receiver};
use rayon::prelude::*;
use std::time::Instant;
use std::time::Duration;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parallel_crossbeam_channel() {
        let (s, r) = bounded(10);
        let mut handles = Vec::new();
        
        for i in 0..10 {
            let s = s.clone();
            let handle = thread::spawn(move || {
                thread::sleep(Duration::from_millis(1));
                s.send(i).unwrap();
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
        drop(s);
        for i in r {
            println!("{}", i);
        }
    }

    #[test]
    fn test_rayon() {
        let mut v = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        v.par_iter_mut().for_each(|x| {
            println!("{}: {:?}", thread::current().name().unwrap_or("2"), x);
            *x += 1
        });
        println!("{:?}", v);
    }

    #[test]
    fn test_complementary() {
        let mut dna_sequence = "AATCGGAATAGAGAGCGATATTACAGAATGAGTAGA".chars().collect::<Vec<char>>();
        dna_sequence.reverse();
        dna_sequence.par_iter_mut().for_each(|x| {
            match x {
                'A' => *x = 'T',
                'T' => *x = 'A',
                'C' => *x = 'G',
                'G' => *x = 'C',
                _ => (),
            }
        });

        let complementary_sequence = dna_sequence.iter().collect::<String>();
        println!("{}", complementary_sequence);
    }
}