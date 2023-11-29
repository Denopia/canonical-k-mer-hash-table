//#include <iostream>
#include <fstream>
#include <cstdint>
#include "functions_strings.hpp"

#pragma once

/*
    Class for reading fasta files. 
    Absolute jank... 
    but should work if the fasta file is properly formatted.
*/

class FastaReader
{
    public:
        std::ifstream fasta_file;
        std::string fasta_path;

        std::string current_read;
        std::string current_read_reversed;
        
        std::string current_line;
        
        uint64_t current_line_number;
        uint64_t current_read_number;
        uint64_t current_read_length;
        bool i_have_reads;
        bool reverse_reads_enabled;
        bool current_is_forward;
        
        // Constructor
        FastaReader(std::string path, bool reverse);

        // Makes the reader roll on to the next read in the input file
        void roll_to_next_read();

        // Reads lines from the input file and builds a single string that contains the current full read
        void read_the_next_read();

        int get_current_read_length();

        bool read_is_loaded();

        char get_current_read_character_at(int position);

        std::string get_current_read();
};