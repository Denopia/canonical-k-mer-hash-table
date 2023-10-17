# Kaarme

Kaarme is a memory-efficient hash table implementation for long k-mers. In this program its functionality is demonstrated with a k-mer counter.

Kaarme k-mer counter is (partially) multithreaded, and counts only canonical k-mers.

Supported input types are fasta and plain text (one read per line) files.

## Requirements
* CMake 3.10
* C++20

(This program has been tested only on Linux based OS)

## Installation

Run the following commands:

```
git clone git@github.com:Denopia/kaarme.git
cd kaarme
mkdir build
cd build
cmake -S .. -B .
cmake --build .
```

## Usage

Run the following command:

```
./kaarme [parameters]
```

Parameters are:
```
-m : Hash table type as integer. Use 0 for plain hash table and 2 for Kaarme hash table.
-i : Input type as integer. Use 0 for fasta and 2 for plain text.
-k : k-mer length as integer.
-s: Hash table size as integer. (Resizing is not implemented at the moment so if the hash table is too small, the program must be restarted manually with bigger hash table size.)
-b: Estimated number of unique k-mers. (Used as the bloom filter size)
-f: Bloom filter false positive rate as double. 
-t: Number of threads as integer. Minimum number of threads is 3.
-a: Minimum numner of k-mer occurrences for it to printed in the output file.
-p: Path to the input file as string.  
-o: Path to the output file as string.
```

Kaarme hash table has two modes: one that uses bloom filter to try to filter out k-mers that occur less than twice, and other that does not use a bloom filter. In the case you want to use bloom filter, you must provide an estimate for the number of unique k-mers in the data set as parameter -b. If you do not wish to use the bloom filter, you must provide a size for the hash table as parameter -s. If -b is provided bloom filter mode is used, otherwise the non-bloom filter mode is used. If neither parameter is provided, the program fails to run correctly. If bloom filter mode is used, you can provide the false positive rate as parameter -f. Otherwise, default value 0.01 is used.

All other parameters are required.

## Example

Use the installation instructions to install the program. Then run the following (assuming you are in the project root directory and Kaaarme is installed in build directory):
```
./build/kaarme -m 2 -i 0 -k 51 -s 8000000 -t 5 -a 2 -p example/ecoli1x.fasta -o example/ecoli1x-51mers.txt
```
Now the example directory should contain a file called ecoli1x-51mers.txt with all 51-mers that appear at least twice in ecoli1x.fasta.

The example above uses the Kaarme without bloom filtering. To use bloom filtering, run the following example:
```
./build/kaarme -m 2 -i 0 -k 51 -b 4000000 -t 5 -a 2 -p example/ecoli1x.fasta -o example/ecoli1x-51mers.txt
```

## Licence

TBD
