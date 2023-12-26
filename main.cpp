#include <iostream>
#include "parallel_parser.hpp"
#include "hash_functions.hpp"
#include "double_bloomfilter.hpp"

#include <cmath>
#include <chrono>

#include "external/CLI11.hpp"
#include <filesystem>
#include <zlib.h>

/*
###########################################################################################################################
This is all I want to see in main(?)
###########################################################################################################################
*/

bool is_gz(std::string& in_file){
    std::ifstream ifs(in_file, std::ios::binary | std::ios::in);
    uint8_t byte1, byte2;
    ifs.read((char *)&byte1, 1);
    ifs.read((char *)&byte2, 2);
    return (byte1 == 0x1f) && (byte2 == 0x8b);
}

std::tuple<char, bool, bool> file_format(std::string& in_file){
    std::tuple<char, bool, bool> res;
    std::get<2>(res) = is_gz(in_file);

    std::filesystem::path pt = std::filesystem::path(in_file);
    std::string ext = pt.extension();
    char sym;

    if(std::get<2>(res)){//gzipped
        while(ext==".gz"){
            pt.replace_extension();
            ext = pt.extension();
        }
        gzFile g_file;
        g_file = gzopen(in_file.c_str(), "r");
        gzread(g_file, &sym, 1);
    } else {
        std::ifstream ifs(in_file, std::ios::binary | std::ios::in);
        ifs.read(&sym, 1);
    }

    bool ill_formed=false;
    if(ext==".fasta" || ext==".fa"){
        if(sym!='>'){
            ill_formed = true;
        }
        std::get<0>(res) = '>';
    }else if(ext==".fastq" || ext==".fq"){
        if(sym!='@'){
            ill_formed = true;
        }
        std::get<0>(res) = '@';
    }else{
        std::string dna_syms = "actgACGT";
        if(std::find(dna_syms.begin(), dna_syms.end(), sym)==dna_syms.end()){
            ill_formed = true;
        }
        std::get<0>(res) = 0;
    }
    std::get<1>(res) = ill_formed;
    return res;
}

struct arguments{
    int hash_table_mode = -1;
    int input_mode = -1;
    uint8_t header_symbol = 0;
    off_t k = 0;
    uint64_t min_slots = 0;
    uint64_t min_abundance = 0;
    size_t n_threads = 1;//number of threads

    std::string input_file;
    std::string output_file;

    bool verbose = false;
    bool user_wants_help = false;
    bool debug = false;

    uint64_t bloom_filter_1_size = 1000000000;
    uint64_t bloom_filter_2_size = 1000000000;
    uint64_t bf1hfn = 1;
    uint64_t bf2hfn = 1;
    double fpr = 0.01;
    uint64_t expected_number_of_unique_kmers = 0;
    bool use_bloom_filter = false;

    bool ver{};
    std::string version ="0.0.1v";
};


//class MyFormatter : public CLI::Formatter {
//public:
//    MyFormatter() : Formatter() {}
//    std::string make_option_opts(const CLI::Option *) const override { return ""; }
//};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option * opt) const override {
        std::stringstream out;
        if(!opt->get_option_text().empty()) {
            out << " " << opt->get_option_text();
        } else {
            if(opt->get_type_size() != 0) {
                if(!opt->get_type_name().empty()) out << " " << get_label(opt->get_type_name());
                //if(!opt->get_default_str().empty()) out << "=" << opt->get_default_str();
                if(opt->get_expected_max() == CLI::detail::expected_max_vector_size) out << " ...";
                else if(opt->get_expected_min() > 1) out << " x " << opt->get_expected();

                if(opt->get_required())
                    out << " " << get_label("REQUIRED");
            }
        }
        return out.str();
    }
};

int parse_app(CLI::App& app, struct arguments& args) {

    auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(29);
    app.formatter(fmt);

    app.add_option("INPUT", args.input_file, "Input file (automatic format detection)")->check(CLI::ExistingFile)->required();
    app.add_option("KLEN", args.k, "k-mer length")->check(CLI::PositiveNumber)->required();

    app.add_option("-m,--hash-table-type", args.hash_table_mode, "Hash table type: 0 for plain and 2 for kaarme (def. 2)")->check(CLI::Range(0,2))->default_val(2);
    app.add_option("-a,--min-k-abu", args.min_abundance, "Minimum abundance threshold for the output k-mers (def. 2)")->default_val(2);
    app.add_option("-t,--threads", args.n_threads, "Number of working threads (def. 3)")->check(CLI::Range(3,64));
    app.add_option("-o,--output-file", args.output_file, "Output file where the k-mer counts will be stored");
    auto *bf_flag = app.add_flag("-b,--use-bfilter", args.use_bloom_filter, "Use bloom filters to discard unique k-mers");
    auto fpr = app.add_option("-f,--bfilter-fpr", args.fpr, "Bloom filter false positive rate (def. 0.01)")->check(CLI::Range(0.001,0.999))->default_val(0.01);

    auto ex_group = app.add_option_group("dummy group2");
    auto *ht_size = ex_group->add_option("-s,--hash-tab-size", args.min_slots, "Hash table size");
    auto bf_unq_kmers = ex_group->add_option("-u,--unq-kmers", args.expected_number_of_unique_kmers, "Estimated number of unique k-mers");
    ex_group->require_option(1);

    bf_flag->needs(bf_unq_kmers);
    bf_unq_kmers->needs(bf_flag);
    fpr->needs(bf_flag);

    ht_size->group("Mandatory params");
    bf_unq_kmers->group("Mandatory params");
    return 0;
}

int main(int argc, char const* argv[])
{
    arguments args;
    CLI::App app("Space-efficient k-mer counter");

    parse_app(app, args);
    CLI11_PARSE(app, argc, argv);

    auto format = file_format(args.input_file);

    if(std::get<1>(format)){
        std::cerr<<"Input file "<<args.input_file<<" is ill-formed"<<std::endl;
        exit(1);
    }

    args.header_symbol = std::get<0>(format);
    bool is_gzipped = std::get<2>(format);

    std::string file = std::filesystem::path(args.input_file).filename();
    std::string fmt;
    if(args.header_symbol=='>'){
        fmt="FASTA";
        args.input_mode = 0;
    }else if(args.header_symbol=='@'){
        fmt="FASTQ";
        args.input_mode = 1;
    }else{
        fmt="ONE-STR-PER-LINE";
        args.input_mode = 2;
    }
    if(args.output_file.empty()){
        args.output_file = std::filesystem::path(args.input_file).replace_extension().filename().string()+".karme_counts";
    }

    std::cout<<"Running settings: "<<std::endl;
    std::cout<<"  input file:               "<<file<<std::endl;
    std::cout<<"  input format:             "<<fmt<<std::endl;
    std::cout<<"  gzip compressed:          "<<(is_gzipped?"yes":"no")<<std::endl;
    std::cout<<"  k-mer length:             "<<args.k<<std::endl;
    std::cout<<"  min. abundance threshold: "<<args.min_abundance<<std::endl;
    std::cout<<"  hash table type:          "<<(args.hash_table_mode==0?"plain":"kaarme")<<std::endl;
    std::cout<<"  using bloom filers:       "<<(args.use_bloom_filter?"yes":"no")<<std::endl;
    if(args.use_bloom_filter){
        std::cout<<"    est. unique k-mers:     "<<args.expected_number_of_unique_kmers<<std::endl;
        std::cout<<"    false positive rate:    "<<args.fpr<<std::endl;
    }else{
        std::cout<<"    est. hash table size:   "<<args.min_slots<<std::endl;
    }
    std::cout<<"  working threads:          "<<args.n_threads<<std::endl;
    std::cout<<"  output file:              "<<args.output_file<<std::endl;

    //if(args.ver){
    //    std::cout<<args.version<<std::endl;
    //    exit(0);
    //}

    /*int hash_table_mode = -1;
    int input_mode = -1;
    uint8_t header_symbol = 0;
    off_t k = 0;
    uint64_t min_slots = 0;
    uint64_t min_abundance = 0;
    size_t n_threads = 3;//number of threads

    std::string input_file;
    std::string output_file;

    bool verbose = false;
    bool user_wants_help = false;
    bool debug = false;

    uint64_t bloom_filter_1_size = 1000000000;
    uint64_t bloom_filter_2_size = 1000000000;
    uint64_t bf1hfn = 1;
    uint64_t bf2hfn = 1;
    double fpr = 0.01;

    uint64_t expected_number_of_unique_kmers = 0;
    bool use_bloom_filter = false;

    int argi = 1;
    // --- Parse arguments ---
    while (argi < argc)
    {
        std::string as(argv[argi]);
        
        if (as.compare("-m") == 0){
            hash_table_mode = stoi(std::string(argv[argi+1]));
            argi += 2;
        }
        else if (as.compare("-b") == 0)
        {
            expected_number_of_unique_kmers = stoi(std::string(argv[argi+1]));
            use_bloom_filter = true;
            argi += 2;
        }
        else if (as.compare("-f") == 0)
        {
            fpr = stod(std::string(argv[argi+1]));
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
        std::cout << "\t-b [integer] bloom filter size";
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

    // Confirm parameters are acceptable
    if ((hash_table_mode < 0) || (hash_table_mode > 2))
    {
        std::cout << "Hash table mode not supported\n";
        exit(1);
    }

    if (args.input_mode == 0){
        args.header_symbol = '>';
    } else if (args.input_mode == 1){
        args.header_symbol = '@';
    } else if (args.input_mode == 2){
        args.header_symbol = 0;
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
        //std::cout << "You must set the size of the hash table\n";
        //exit(1);
    }

    if (n_threads < 3)
    {
        std::cout << "Number of threads must be at least 3\n";
        exit(1);
    }

    //if (min_abundance < 0)
    //{
    //    std::cout << "Minimum k-mer frequency for outputting must be set\n";
    //    exit(1);
    //}

    if (input_file.empty())
    {
        std::cout << "Path to input file was not given\n";
        exit(1);
    }

    if (output_file.empty())
    {
        std::cout << "Path to output file was not given\n";
        exit(1);
    }*/

    args.n_threads = args.n_threads - 2;

    //settings for the file buffers
    size_t active_chunks = 2*args.n_threads; //number of chunks in the buffer
    off_t chunk_size = 1024*1024*10; //size in bytes for every chunk

    //=============================================================================================================================================
    //=============================================================================================================================================
    //        Do bloom filtering here
    //=============================================================================================================================================
    //=============================================================================================================================================
    
    if(args.use_bloom_filter)
    {
        // Set bloom filter threads
        uint64_t bf_threads = args.n_threads;
        //uint64_t bf_threads = 1;

        double bloom_filter_error_rate = args.fpr;
        double bloom_filter_bits_min = (-double(args.expected_number_of_unique_kmers) * std::log(bloom_filter_error_rate)) / (std::pow(std::log(2), 2));
        double hash_functions = (bloom_filter_bits_min / double(args.expected_number_of_unique_kmers)) * std::log(2);
        uint64_t bloom_filter_bits_2P = 2;

        while (bloom_filter_bits_2P < uint64_t(bloom_filter_bits_min)){
            bloom_filter_bits_2P = 2*bloom_filter_bits_2P;
        }

        uint64_t bloom_filter_bits = bloom_filter_bits_2P;
        //uint64_t bloom_filter_bits = bloom_filter_bits_min;

        // Bloom filter sizes
        args.bloom_filter_1_size = bloom_filter_bits;
        //bloom_filter_2_size = bloom_filter_bits;
        // Number of Bloom filter hash functions
        args.bf1hfn = std::ceil(hash_functions);
        args.bf2hfn = std::ceil(hash_functions);

#ifdef DEBUG
        std::cout << "Bloom filter error rate " << bloom_filter_error_rate << "\n";
        std::cout << "Bloom filter bits " << args.bloom_filter_1_size << "\n";
        std::cout << "Number of hash functions " << args.bf1hfn << "\n";
#endif
    
        // Set sizes
        uint64_t bf1_size = args.bloom_filter_1_size;
        //uint64_t bf2_size = bloom_filter_1_size;

        //const uint64_t hffs = 7;
        auto * double_adbf = new DoubleAtomicDoubleBloomFilter(bf1_size, args.bf1hfn);
       
        uint64_t rolling_hasher_mod = uint64_t(1) << 54;
        uint64_t bf1_multiplier = 5;
        uint64_t bf1_modmulinv = mathfunctions::modular_multiplicative_inverse_coprimes(bf1_multiplier, rolling_hasher_mod);

#ifdef DEBUG
        std::cout << "Modular multiplicative inverse for A=" << bf1_multiplier << " and M=" << rolling_hasher_mod << " is " << bf1_modmulinv << "\n";
#endif

        // Bloom filter Atomic Double Bloom Filter
        parse_input_pointer_atomic_variable_BLOOM_FILTERING<uint8_t, false>()(bf1_modmulinv, bf1_multiplier, double_adbf, bf1_size,
                                                                    rolling_hasher_mod, args.bf1hfn,
                                                                    args.input_file, chunk_size, active_chunks, bf_threads, args.k,
                                                                    args.header_symbol, args.input_mode, args.debug);

        //auto end_bf = std::chrono::high_resolution_clock::now();
        //auto duration_bf = std::chrono::duration_cast<std::chrono::microseconds>(end_bf - start_bf);
        //std::cout << "Time used for Bloom filtering: " << duration_bf.count() << " microseconds\n";

        //exit(0);

        //if(args.use_bloom_filter){
        args.min_slots = 2*double_adbf->get_new_in_second();
        //}

#ifdef DEBUG
        std::cout << "... Resizing bloom filter ...\n";
#endif
        auto start_resizing = std::chrono::high_resolution_clock::now();
        double_adbf->resize();
        auto end_resizing = std::chrono::high_resolution_clock::now();
        auto resizing_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_resizing - start_resizing);
#ifdef DEBUG
        std::cout << "Time used to resize bloom filter: " << resizing_duration.count() << " microseconds\n";
#endif

        if(args.hash_table_mode == 0)
        {
            if(is_gzipped){
                parse_input_atomic_flag_BF<uint8_t, true>()(bf1_modmulinv, bf1_multiplier, double_adbf, bf1_size, 
                                                                    rolling_hasher_mod, hash_functions,
                                                                    args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k,
                                                                    args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }else{
                parse_input_atomic_flag_BF<uint8_t, false>()(bf1_modmulinv, bf1_multiplier, double_adbf, bf1_size, 
                                                                    rolling_hasher_mod, hash_functions,
                                                                    args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k,
                                                                    args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }    
        }
        else if (args.hash_table_mode == 1)
        {
            if(is_gzipped){
                parse_input_pointer_atomic_flag<uint8_t, true>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }else{
                parse_input_pointer_atomic_flag<uint8_t, false>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }
        }
        else if (args.hash_table_mode == 2)
        {
            if(is_gzipped){
                parse_input_pointer_atomic_variable_BF<uint8_t, true>()(bf1_modmulinv, bf1_multiplier, double_adbf, bf1_size, 
                                                                    rolling_hasher_mod, hash_functions, 
                                                                    args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k,
                                                                    args.header_symbol, args.min_slots, args.min_abundance, args.input_mode, args.debug);
            }else{
                parse_input_pointer_atomic_variable_BF<uint8_t, false>()(bf1_modmulinv, bf1_multiplier, double_adbf, bf1_size, 
                                                                    rolling_hasher_mod, hash_functions,
                                                                    args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k,
                                                                    args.header_symbol, args.min_slots, args.min_abundance, args.input_mode, args.debug);
            }
        }
        else
        {
            std::cout << "Chosen mode not recognized\n";
        }

        delete double_adbf;
    }

    // If bloom filter is not used
    else
    {
        if (args.hash_table_mode == 0)
        {
            if(is_gzipped){
                parse_input_atomic_flag<uint8_t, true>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }else{
                parse_input_atomic_flag<uint8_t, false>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }    
        }
        else if (args.hash_table_mode == 1)
        {
            if(is_gzipped){
                parse_input_pointer_atomic_flag<uint8_t, true>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }else{
                parse_input_pointer_atomic_flag<uint8_t, false>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode);
            }
        }
        else if (args.hash_table_mode == 2)
        {
            if(is_gzipped){
                parse_input_pointer_atomic_variable<uint8_t, true>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode, args.debug);
            }else{
                parse_input_pointer_atomic_variable<uint8_t, false>()(args.input_file, args.output_file, chunk_size, active_chunks, args.n_threads, args.k, args.header_symbol, args.min_slots, args.min_abundance, args.input_mode, args.debug);
            }
        }
        else
        {
            std::cout << "Chosen mode not recognized\n";
        }
    }
}