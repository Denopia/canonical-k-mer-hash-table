# Kaarme

Kaarme is a memory-efficient hash table implementation for long k-mers. In this program its functionality is demonstrated with a k-mer counter.

Kaarme k-mer counter is (partially) multithreaded, and counts only canonical k-mers.

Supported input types are fasta and plain text (one read per line) files.

## Requirements
* CMake 3.10
* C++17

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
./kaarme [OPTIONS] INPUT KLEN TABLE_SIZE
```

Parameters are:
```

Positionals:
  INPUT TEXT REQUIRED        Input file
  KLEN INT REQUIRED          k-mer length
  TABLE_SIZE UINT REQUIRED   Hash table size

Options:
  -h,--help                  Print this help message and exit
  -m,--hash-table-type INT   Hash table type: 0 for plain and 2 for kaarme (def. 2)
  -a,--min-k-abu UINT        Minimum abundance threshold for the output k-mers (def. 2)
  -t,--threads UINT          Number of hashing threads (def. 1)
  -o,--output-file TEXT      Output file where the k-mer counts will be stored

Bloom filter options:
  -b,--use-bfilter           Use bloom filters to discard low-frequency k-mers
  -u,--unq-kmers UINT        Estimated number of unique k-mers
  -f,--bfilter-fpr FLOAT     Bloom filter false positive rate (def. 0.01)
```

Kaarme hash table has two modes: one that uses Bloom filter to try to filter out k-mers that occur less than twice, and other that does not use a Bloom filter. In the case you want to use Bloom filter, you must provide an estimate for the number of unique k-mers in the data set as parameter -b. If you do not wish to use the Bloom filter, you must provide a size for the hash table as parameter -s. If -b is provided Bloom filter mode is used, otherwise the non-Bloom filter mode is used. If neither parameter is provided, the program fails to run correctly. If Bloom filter mode is used, you can provide the false positive rate as parameter -f. Otherwise, default value 0.01 is used.

All other parameters are required.

## Example

Use the installation instructions to install the program. Then run the following (assuming you are in the project root directory and Kaaarme is installed in build directory):
```
./build/kaarme example/ecoli1x.fasta 51 8000000 -t 3
```
Now the example directory should contain a file called ecoli1x-51mers.txt with all 51-mers that appear at least twice in ecoli1x.fasta.

The example above uses the Kaarme without Bloom filtering. To use Bloom filtering, run the following example:
```
./build/kaarme example/ecoli1x.fasta 51 8000000 -t 3 -u 4000000 --use-bfilter 
```

This implementation also includes a basic k-mer counter using a plain hash table instead of the Kaarme hash table. This plain hash table k-mer counter can be used by setting the hash table type parameter -m value to 0. To run the above examples using the plain hash table, run the following without Bloom filter:
```
./build/kaarme example/ecoli1x.fasta 51 8000000 -t 3 -m 0
```
And the following for plain hash table with Bloom filter:
```
./build/kaarme example/ecoli1x.fasta 51 8000000 -t 3 -m 0 -u 4000000 --use-bfilter
```

## Licence

TBD
