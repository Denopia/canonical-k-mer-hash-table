#include <iostream>
//#include "program_runs.hpp"
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
    if(argc!=8){
        std::cout<<"Some error"<<std::endl;
        exit(0);
    }

    //bool use_atomic_flag = true;

    int mode = -1;

    //settings for the file buffers
    size_t n_threads= stoi(std::string(argv[4]));//number of threads
    size_t active_chunks = 10; //number of chunks in the buffer
    off_t chunk_size = 1024*1024*10; //size in bytes for every chunk

    std::string input_file = std::string(argv[7]);
    std::string output_file = std::string(argv[6]);
    off_t k=25;//a value to test the kmers
    bool is_gzipped = false;//TODO check if its gzipped

    k = stoi(std::string(argv[2]));
    uint64_t min_slots = stoi(std::string(argv[3]));
    uint64_t min_abundance = stoi(std::string(argv[5]));


    // Choose mode
    // 0 for basic hash table with atomic flags
    if (stoi(std::string(argv[1])) == 0)
        mode = 0;
    // 1 for pointer hash table with atomic variables
    else if (stoi(std::string(argv[1])) == 1)
        mode = 1;
    else if (stoi(std::string(argv[1])) == 2)
        mode = 2;
    else
        mode = -1;

    

    if (mode == 0)
    {
        if(is_gzipped){
            parse_input_atomic_flag<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }else{
            parse_input_atomic_flag<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }    
    }
    else if (mode == 1)
    {
        if(is_gzipped){
            parse_input_pointer_atomic_flag<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }else{
            parse_input_pointer_atomic_flag<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }
    }
    else if (mode == 2)
    {
        if(is_gzipped){
            parse_input_pointer_atomic_variable<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }else{
            parse_input_pointer_atomic_variable<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, min_slots, min_abundance);
        }
    }
    else
    {
        std::cout << "Chosen mode not recognized\n";
    }
    
}
