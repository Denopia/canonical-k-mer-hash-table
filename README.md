# Kaarme

Kaarme is a memory-efficient hash table implementation for long k-mers. In this software, its functionality is demonstrated with a k-mer counter.

Kaarme k-mer counter is (partially) multithreaded, and counts only canonical k-mers.

Supported input types are fasta and plain text (one read per line) files.

## Requirements
* CMake 3.10
* C++17

This program has been tested only on OSx and Linux-based OS

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
./kaarme [OPTIONS] INPUT KLEN TABLE_SIZE
```

Parameters are:
```
Positionals:
  INPUT TEXT REQUIRED        Input file
  KLEN INT REQUIRED          k-mer length

Options:
  -h,--help                  Print this help message and exit
  -m,--hash-table-type INT   Hash table type: 0 for plain and 2 for kaarme (def. 2)
  -a,--min-k-abu UINT        Minimum abundance threshold for the output k-mers (def. 2)
  -t,--threads UINT          Number of working threads (def. 3)
  -o,--output-file TEXT      Output file where the k-mer counts will be stored
  -b,--use-bfilter           Use bloom filters to discard unique k-mers
  -f,--bfilter-fpr FLOAT     Bloom filter false positive rate (def. 0.01)


[Exactly 1 of the following options is required]
Excluding params:
  -s,--hash-tab-size UINT    Hash table size
  -u,--unq-kmers UINT        Estimated number of unique k-mers
```

Kaarme hash table has two modes: one that uses Bloom filter to filter out k-mers that occur less than twice
(--use-bfilter flag), and other that does not use a Bloom filter. In the case you want to use Bloom filter, you must
provide an estimate for the number of unique k-mers (--unq-kmers parameter). If you do not wish to use the Bloom
filter, you must provide a size for the hash table as parameter -s. If Bloom filter mode is used, you can provide the
false positive rate as parameter -f. Otherwise, default value 0.01 is used.

## Example

Use the installation instructions to install the program. Then run the following (assuming you are in the project root directory and Kaarme is installed in build directory):
```
./build/kaarme example/ecoli1x.fasta 51 -s 8000000 -t 3 -o example/ecoli1x-51mers.txt
```
Now the example directory should contain a file called ecoli1x-51mers.txt with all 51-mers that appear at least twice in ecoli1x.fasta.

The example above uses the Kaarme without Bloom filtering. To use Bloom filtering, run the following example:
```
./build/kaarme example/ecoli1x.fasta 51 -t 3 -u 4000000 --use-bfilter -o example/ecoli1x-51mers.txt 
```

This implementation also includes a basic k-mer counter using a plain hash table instead of the Kaarme hash table. This plain hash table k-mer counter can be used by setting the hash table type parameter -m value to 0. To run the above examples using the plain hash table, run the following without Bloom filter:
```
./build/kaarme example/ecoli1x.fasta 51 -s 8000000 -t 3 -m 0 -o example/ecoli1x-51mers.txt
```
And the following for plain hash table with Bloom filter:
```
./build/kaarme example/ecoli1x.fasta 51 -t 3 -m 0 -u 4000000 --use-bfilter -o example/ecoli1x-51mers.txt
```

## Licence

TBD
