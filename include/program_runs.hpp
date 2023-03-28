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


// Mode for canonical hash table
int run_mode_1(int argc, char const* argv[]);
