# Kaarme

Kaarme is a memory-efficient hash table implementation for long k-mers. In this program its functionality is demonstrated with a k-mer counter.

Kaarme k-mer counter is (partially) multithreaded, and counts only canonical k-mers.

Supported input types are fasta and plain text (one read per line) files.

## Requirements
CMake 3.10

C++20

## Installation

Run the following commands (in the directory where this README file is):

1. mkdir build
2. cd build
3. cmake -S ../source -B .
4. cmake --build .

## Usage

Run the following command:

./kht [parameters]

Parameters are:
- -m : Hash table type as integer. Use 0 for plain hash table and 2 for Kaarme hash table.
- -i : Input type as integer. Use 0 for fasta and 2 for plain text.
- -k : k-mer length as integer.
- -s: Hash table size as integer. (Resizing is not implemented at the moment so if the hash table is too small, the program must be restarted manually with bigger hash table size.)
- -t: Number of threads as integer. Minimum number of threads is 3.
- -a: Minimum numner of k-mer occurrences for it to printed in the output file.
- -p: Path to the input file as string.  
- -o: Path to the output file as string.

All parameters are required.

## Licence

TBD
