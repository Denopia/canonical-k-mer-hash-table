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
    //if(argc!=9){
    //    std::cout<<"Some error"<<std::endl;
    //    exit(0);
    //}

    int hash_table_mode = -1;
    int input_mode = -1;
    uint8_t header_symbol = 0;
    off_t k = 0;
    uint64_t min_slots = 0;
    uint64_t min_abundance = 0;
    size_t n_threads = 3;//number of threads

    std::string input_file = "";
    std::string output_file = "";

    bool verbose = false;
    bool user_wants_help = false;
    bool debug = false;

    int argi = 1;
    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        
        if (as.compare("-m") == 0){
            hash_table_mode = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-i") == 0)
        {
            input_mode = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-k") == 0)
        {
            k = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-s") == 0)
        {
            min_slots = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-t") == 0)
        {
            n_threads = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-a") == 0)
        {
            min_abundance = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-p") == 0)
        {
            input_file = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-o") == 0)
        {
            output_file = std::string(argv[argi+1]);
            argi += 2;
        }
        else if (as.compare("-v") == 0)
        {
            verbose = true;
            argi += 1;
        }
        else if (as.compare("-d") == 0)
        {
            debug = true;
            argi += 1;
        }
        else
        {
            argi++;
        }
    }

    if (user_wants_help)
    {   std::cout << "This is a k-mer has table program. Here are the accepted arguments:\n\n";
        std::cout << "\t-k [integer] determines the k-mer length\n";
        std::cout << "\t-a [integer] minimum k-mer abundance (for writing output, default: 1)\n";
        std::cout << "\t-s [integer] determines the minimum hash table size (actual used size is the next prime)\n";
        std::cout << "\t-u [integer] determines a limit for the number of items that can be stored in the hash table (default: maximum size)\n";
        std::cout << "\t-q flag enables queries after hash table is built\n";
        std::cout << "\t-f flag tells the program to do a first k-mer check after buffer is empty (should reduce the number of stored entries)\n";
        std::cout << "\t-p [string] indicates the path to the read file (at the moment fasta format is required)\n";
        std::cout << "\t-v flag tells the program to print what is going on during the program\n";
        std::cout << "\t-o [string] tells the program to write found k-mers into an output file with this name\n";
        std::cout << "\t-d [string] shows the directory where output is written\n";
        std::cout << "\t-h flag shows this help message\n\n";

        return 0;
    }

    // Confrim parameters are acceptable
    if ((hash_table_mode < 0) || (hash_table_mode > 2))
    {
        std::cout << "Hash table mode not supported\n";
        exit(1);
    }

    if (input_mode == 0){
        header_symbol = '>';
    } else if (input_mode == 1){
        header_symbol = '@';
    } else if (input_mode == 2){
        header_symbol = 0;
    } else {
        std::cout << "Input mode not supported\n";
        exit(1);
    }

    if (k == 0)
    {
        std::cout << "Value of k must be set and greater than 0\n";
        exit(1);
    }

    if (min_slots == 0)
    {
        std::cout << "You must set the size of the hash table\n";
        exit(1);
    }

    if (n_threads < 3)
    {
        std::cout << "Number of threads must be at least 3\n";
        exit(1);
    }

    if (min_abundance == 0)
    {
        std::cout << "Minimum k-mer frequency for outputting must be set\n";
        exit(1);
    }

    if (input_file == "")
    {
        std::cout << "Path to input file was not given\n";
        exit(1);
    }

    if (output_file == "")
    {
        std::cout << "Path to output file was not given\n";
        exit(1);
    }

    //settings for the file buffers
    size_t active_chunks = 12; //number of chunks in the buffer
    off_t chunk_size = 1024*1024*10; //size in bytes for every chunk

    bool is_gzipped = false;//TODO check if its gzipped

    if (hash_table_mode == 0)
    {
        if(is_gzipped){
            parse_input_atomic_flag<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode);
        }else{
            parse_input_atomic_flag<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode);
        }    
    }
    else if (hash_table_mode == 1)
    {
        if(is_gzipped){
            parse_input_pointer_atomic_flag<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode);
        }else{
            parse_input_pointer_atomic_flag<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode);
        }
    }
    else if (hash_table_mode == 2)
    {
        if(is_gzipped){
            parse_input_pointer_atomic_variable<uint8_t, true>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode, debug);
        }else{
            parse_input_pointer_atomic_variable<uint8_t, false>()(input_file, output_file, chunk_size, active_chunks, n_threads, k, header_symbol, min_slots, min_abundance, input_mode, debug);
        }
    }
    else
    {
        std::cout << "Chosen mode not recognized\n";
    }
    
}
