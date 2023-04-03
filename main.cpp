#include <iostream>
#include "program_runs.hpp"
#include "parallel_parser.hpp"

/*
###########################################################################################################################
This is all I want to see in main(?)
###########################################################################################################################
*/

int main(int argc, char const* argv[])
{
    // Canonical hash table mode
    //run_mode_1(argc, argv);
    if(argc!=2){
        std::cout<<"Some error"<<std::endl;
        exit(0);
    }

    //settings for the file buffers
    size_t n_threads= 4;//number of threads
    size_t active_chunks = 10; //number of chunks in the buffer
    off_t chunk_size = 1024*1024*20; //size in bytes for every chunk

    std::string input_file = std::string(argv[1]);
    off_t k=25;//a value to test the kmers
    bool is_gzipped = false;//TODO check if its gzipped

    if(is_gzipped){
        parse_input<uint8_t, true>()(input_file, chunk_size, active_chunks, n_threads, k);
    }else{
        parse_input<uint8_t, false>()(input_file, chunk_size, active_chunks, n_threads, k);
    }
}
