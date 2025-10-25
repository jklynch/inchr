use anyhow::{bail, Result};
use clap::Parser;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::io::Write;
use std::time::{Duration, Instant};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// FASTQ file path
    #[arg(short, long)]
    fastq: String,

    /// k-mer length
    #[arg(short, long)]
    kmer_length: usize,

    /// Number of top k-mers to show
    #[arg(short, long, default_value_t = 10)]
    top: usize,

    /// Output FASTQ file path
    #[arg(short, long)]
    output: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct KmerInfo {
    pub count: u64,
    pub entropy: f64,
}


fn find_and_count_kmers(
    fastq_path: &str,
    kmer_length: usize,
) -> Result<(HashMap<Vec<u8>, KmerInfo>, Duration)> {
    let mut reader = parse_fastx_file(fastq_path)?;
    let mut seqs = Vec::new();
    while let Some(record) = reader.next() {
        seqs.push(record?.seq().to_vec());
    }

    let start_time = Instant::now();

    let mut kmer_counts: HashMap<Vec<u8>, KmerInfo> = HashMap::new();

    for seq in seqs {
        for kmer in seq.windows(kmer_length) {
            let entry = kmer_counts.entry(kmer.to_vec()).or_insert(KmerInfo { count: 0, entropy: 0.0 });
            entry.count += 1;

            // Calculate entropy for the k-mer
            let mut nucleotide_counts = HashMap::new();
            for &base in kmer {
                *nucleotide_counts.entry(base).or_insert(0) += 1;
            }

            let mut entropy = 0.0;
            for (_, &count) in &nucleotide_counts {
                let probability = count as f64 / kmer_length as f64;
                if probability > 0.0 {
                    entropy -= probability * probability.log2();
                }
            }
            entry.entropy = entropy;
        }
    }

    let elapsed_time = start_time.elapsed();
    Ok((kmer_counts, elapsed_time))
}

fn inchworm_assemble(kmer_table: &mut HashMap<Vec<u8>, KmerInfo>, kmer_length: usize) -> Vec<u8> {
    if kmer_table.is_empty() {
        return Vec::new();
    }

    // 1) select the k-mer with the highest count to be the seed
    let (seed, _) = kmer_table
        .iter()
        .max_by_key(|&(_, kmer_info)| kmer_info.count)
        .unwrap();
    let mut assembled_sequence = seed.clone();
    kmer_table.remove(&assembled_sequence);

    // Right extension
    loop {
        let mut best_candidate = None;
        let mut max_count = 0;

        for nucleotide in [b'A', b'C', b'G', b'T'].iter() {
            let mut candidate = assembled_sequence.clone();
            candidate.push(*nucleotide);
            let right_most_kmer = &candidate[candidate.len() - kmer_length..];

            if let Some(kmer_info) = kmer_table.get(right_most_kmer) {
                if kmer_info.count > max_count {
                    max_count = kmer_info.count;
                    best_candidate = Some(candidate.clone());
                }
            }
        }

        if let Some(best) = best_candidate {
            assembled_sequence = best;
            let right_most_kmer = &assembled_sequence[assembled_sequence.len() - kmer_length..];
            kmer_table.remove(right_most_kmer);
        } else {
            break;
        }
    }

    // Left extension
    loop {
        let mut best_candidate = None;
        let mut max_count = 0;

        for nucleotide in [b'A', b'C', b'G', b'T'].iter() {
            let mut candidate = vec![*nucleotide];
            candidate.extend_from_slice(&assembled_sequence);
            let left_most_kmer = &candidate[..kmer_length];

            if let Some(kmer_info) = kmer_table.get(left_most_kmer) {
                if kmer_info.count > max_count {
                    max_count = kmer_info.count;
                    best_candidate = Some(candidate.clone());
                }
            }
        }

        if let Some(best) = best_candidate {
            assembled_sequence = best;
            let left_most_kmer = &assembled_sequence[..kmer_length];
            kmer_table.remove(left_most_kmer);
        } else {
            break;
        }
    }

    assembled_sequence
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.kmer_length == 0 {
        bail!("k-mer length must be greater than 0");
    }

    let (kmer_counts, elapsed_time) = find_and_count_kmers(&args.fastq, args.kmer_length)?;

    let total_kmers: u64 = kmer_counts.values().map(|info| info.count).sum();
    println!("Total k-mers: {}", total_kmers);
    println!("Time taken to count k-mers: {:?}", elapsed_time);

    let assembled_sequence = inchworm_assemble(&mut kmer_counts.clone(), args.kmer_length);
    println!("\nAssembled sequence:");
    println!("{}", String::from_utf8_lossy(&assembled_sequence));

    let mut output_file = std::fs::File::create(&args.output)?;
    writeln!(output_file, "@assembled_sequence")?;
    writeln!(output_file, "{}", String::from_utf8_lossy(&assembled_sequence))?;
    writeln!(output_file, "+")?;
    writeln!(output_file, "{}", "F".repeat(assembled_sequence.len()))?;

    let mut sorted_kmers: Vec<_> = kmer_counts.into_iter().collect();
    sorted_kmers.sort_by(|a, b| b.1.count.cmp(&a.1.count));

    println!("\nTop {} k-mers:", args.top);
    for (kmer, kmer_info) in sorted_kmers.iter().take(args.top) {
        println!("{}\t{}\t{:.2}", String::from_utf8_lossy(kmer), kmer_info.count, kmer_info.entropy);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;
    use std::collections::HashMap;

    #[test]
    fn test_find_and_count_kmers() {
        let mut tmpfile = NamedTempFile::new().unwrap();
        writeln!(tmpfile, "@seq1").unwrap();
        writeln!(tmpfile, "ACGTACGT").unwrap();
        writeln!(tmpfile, "+").unwrap();
        writeln!(tmpfile, "FFFFFFFF").unwrap();

        let (kmer_counts, _) =
            find_and_count_kmers(tmpfile.path().to_str().unwrap(), 4).unwrap();

        let mut expected_counts = HashMap::new();
        expected_counts.insert(b"ACGT".to_vec(), KmerInfo { count: 2, entropy: 2.0 });
        expected_counts.insert(b"CGTA".to_vec(), KmerInfo { count: 1, entropy: 2.0 });
        expected_counts.insert(b"GTAC".to_vec(), KmerInfo { count: 1, entropy: 2.0 });
        expected_counts.insert(b"TACG".to_vec(), KmerInfo { count: 1, entropy: 2.0 });

        assert_eq!(kmer_counts, expected_counts);
    }

    #[test]
    fn test_inchworm_assemble() {
        let mut kmer_counts = HashMap::new();
        kmer_counts.insert(b"ACGT".to_vec(), KmerInfo { count: 3, entropy: 0.0 });
        kmer_counts.insert(b"CGTA".to_vec(), KmerInfo { count: 2, entropy: 0.0 });
        kmer_counts.insert(b"GTAC".to_vec(), KmerInfo { count: 1, entropy: 0.0 });
        kmer_counts.insert(b"TACG".to_vec(), KmerInfo { count: 4, entropy: 0.0 });
        let _assembled_sequence = inchworm_assemble(&mut kmer_counts, 4);
    }
}