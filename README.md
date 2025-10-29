# inchr

## build
```sh
$ cargo build -r
```

## usage
```sh
$ target/release/inchr --help
Usage: inchr [OPTIONS] --fastq <FASTQ> --kmer-length <KMER_LENGTH> --output <OUTPUT>

Options:
  -f, --fastq <FASTQ>              FASTQ file path
  -k, --kmer-length <KMER_LENGTH>  k-mer length
  -t, --top <TOP>                  Number of top k-mers to show [default: 10]
  -o, --output <OUTPUT>            Output FASTQ file path
  -h, --help                       Print help
  -V, --version                    Print version
```

## example
```sh
$ target/release/inchr --fastq reads.fq --kmer-length 32 --output reads_assembled.fq
Total k-mers: 1375875
Time taken to count k-mers: 153.53468ms

Top 10 k-mers:
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	599	0.00
GGCCGCACCTTGGTCTGAGTCTGTCGGGAAGC	293	1.90
GGGCCGCACCTTGGTCTGAGTCTGTCGGGAAG	293	1.88
GCCGCACCTTGGTCTGAGTCTGTCGGGAAGCT	292	1.92
GGGGCCGCACCTTGGTCTGAGTCTGTCGGGAA	292	1.88
CCGCACCTTGGTCTGAGTCTGTCGGGAAGCTG	291	1.92
ATAAATACGGGGCCGCACCTTGGTCTGAGTCT	290	1.99
TAAATACGGGGCCGCACCTTGGTCTGAGTCTG	290	1.98
CGGGGCCGCACCTTGGTCTGAGTCTGTCGGGA	290	1.84
ATTCCATAAATACGGGGCCGCACCTTGGTCTG	289	1.99

Time taken to assemble sequences: 184.332015ms
Number of assembled sequences: 12606
Sequence 1  : Length = 6308 , First 32 bases = TTGGCTTGTCTTTAAAAGCCCTAAGAGAGTGA
Sequence 2  : Length = 5661 , First 32 bases = AAGGTCTTGAGGTTGATAAATGGAGACCTGTC
Sequence 3  : Length = 3453 , First 32 bases = GGAAAGGAGAGATCAGAAACGTGGAGACCAAA
Sequence 4  : Length = 3249 , First 32 bases = CCCTGCTCCTGGGTCTCCTACACCTGCTTTTG
Sequence 5  : Length = 3097 , First 32 bases = CGAGGTGACAGAAGCCTCTCCTCCTTGGGAGC
Sequence 6  : Length = 2190 , First 32 bases = AATGGATTGAATCTTATTTTGTACAAACTAAA
Sequence 7  : Length = 2094 , First 32 bases = TTTTCTTGCAATACACAAAAGTTTATTTATTT
Sequence 8  : Length = 1486 , First 32 bases = GTGGAGGTGAGGGGTCCCAGGGGTTTCTTGGG
Sequence 9  : Length = 1323 , First 32 bases = CACTTGAGAGCCCAAAATATTTTATTTAACAA
Sequence 10 : Length = 1278 , First 32 bases = CTGGGGCCAGCAGAGGCGATGCTTATCCATGT
```
