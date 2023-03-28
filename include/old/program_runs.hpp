#include <iostream>
#include <cstdint>
#include <math.h>
#include <chrono>
#include "functions_math.hpp"
#include "file_reader.hpp"
#include "kmer_factory.hpp"
#include "hash_functions.hpp"
#include "kmer_hash_table.hpp"

#pragma once

/*

  -  Here we have program run instructions for all different hashmap versions -

*/


/*
    Normal hashmap mode.

    Runs the program with a normal hashmap with linear(?) probing.
    Still, the k-mers are stored as integers and each character takes 2 bits.

*/
int run_mode_0(int argc, char const* argv[]);


/*
    Single character with pointer to previous k-mer hashmap mode.

    Runs the program with the single character plus pointer scheme.
    All k-mer info in a single integer. Pointer takes another integer

*/
int run_mode_1(int argc, char const* argv[]);

/*
    Single character with pointer to previous k-mer hashmap mode.

    Runs the program with the single character plus pointer scheme.
    All k-mer info in a single integer. Pointer takes another integer
    
    Uses temp storage to avoid storing incomplete k-mers 

*/
int run_mode_2(int argc, char const* argv[]);

