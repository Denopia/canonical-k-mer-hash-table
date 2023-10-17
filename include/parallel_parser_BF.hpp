//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_PARALLEL_PARSER_BF_HPP
#define PARALLEL_PARSING_PARALLEL_PARSER_BF_HPP

#include "text_reader.h"
#include "ts_queue.h"
#include "functions_math.hpp"
#include <thread>
#include <zlib.h>
#include <cstring>
#include <bitset>
#include <chrono>
#include <boost/dynamic_bitset.hpp>
#include "functions_bloom_filter.hpp"
#include "double_bloomfilter.hpp"

enum bio_format{FASTA=0, FASTQ=1, PLAIN=2};

// ==============================================================================================================
// ATOMIC VARIABLES VERSION, MODE=2, WITH BLOOM FILTER
// ==============================================================================================================

template<class sym_type,
         bool is_gzipped=false>
struct parse_input_pointer_atomic_variable_BF{

    

    void operator()(uint64_t bf_modmulinv, uint64_t bf_multiplier, DoubleAtomicDoubleBloomFilter * adbf, uint64_t bloom_filter_size, 
                    uint64_t rolling_hasher_mod, uint64_t hash_functions,
                    std::string& input_file,  off_t chunk_size, size_t active_chunks, size_t n_threads, off_t k,
                    sym_type start_symbol, int input_mode, bool debug){
        
        std::cout << "Starting parallel bloom filtering\n";

        bool print_times = true;
        bool print_other_stuff = true;
        auto start_filtering = std::chrono::high_resolution_clock::now();
        // Create the hash table
        uint64_t kmer_len = k;

        using chunk_type = text_chunk<sym_type>;

        ts_queue<size_t> in_queue;// thread-safe queue that manage the chunks that are ready to be used
        ts_queue<size_t> out_queue; // thread-safe queue that stores the chunks that can be reused for new chunks
        std::vector<chunk_type> text_chunks;
        int fd = open(input_file.c_str(), O_RDONLY);

        // this is for later: to manage compressed inputs
        gzFile gfd;
        if constexpr (is_gzipped){//managed at compilation time
            gfd = gzdopen(fd, "r");
        }

        //get the file size
        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return;

        size_t format; //manage to get the input format
        if (input_mode == 2)
            format = PLAIN;
        else if (input_mode == 0)
            format = FASTA;
        else
        {
            std::cout << "Input file format not supported.";
            return;
        }   

        //lambda function that manages IO operations
        //we feed this function to std::thread
        auto io_worker = [&]() -> void {

            off_t rem_bytes = st.st_size;

#ifdef __linux__
            posix_fadvise(fd, 0, rem_bytes, POSIX_FADV_SEQUENTIAL);//tell the linux kernel we will access the file sequentially so it can use the readahead heuristic more effectively
#endif

            size_t chunk_id=0;
            text_chunks.resize(active_chunks);
            off_t tmp_ck_size;
            bool broken_header=false;


            while(chunk_id<active_chunks && rem_bytes>=k){

                tmp_ck_size = std::min(chunk_size, rem_bytes);
                text_chunks[chunk_id].bytes = tmp_ck_size;
                text_chunks[chunk_id].buffer = (sym_type *)malloc(tmp_ck_size);
                text_chunks[chunk_id].id = chunk_id;
                text_chunks[chunk_id].broken_header = broken_header;
                if constexpr (is_gzipped){
                    rem_bytes = read_chunk_from_gz_file<chunk_type>(gfd, text_chunks[chunk_id], rem_bytes, k-1, broken_header, start_symbol);
                }else{
                    rem_bytes = read_chunk_from_file<chunk_type>(fd, text_chunks[chunk_id], rem_bytes, k-1, broken_header, start_symbol);
                }
                in_queue.push(chunk_id);//as soon as we push, the chunks become visible of the worker threads to consume them
                chunk_id++;
            }

            size_t buff_idx;
            while(rem_bytes>=k){
                out_queue.pop(buff_idx);//it will wait until out_strings contains something
                text_chunks[buff_idx].id = chunk_id++;
                text_chunks[buff_idx].broken_header = broken_header;
                if constexpr (is_gzipped){
                    rem_bytes = read_chunk_from_gz_file<chunk_type>(gfd, text_chunks[buff_idx], rem_bytes, k-1, broken_header, start_symbol);
                }else{
                    rem_bytes = read_chunk_from_file<chunk_type>(fd, text_chunks[buff_idx], rem_bytes, k-1, broken_header, start_symbol);
                }
                in_queue.push(buff_idx);
            }

            //wait for the chunks to be fully processed
            while(!in_queue.empty());

            //remove the unused chunks from the out queue
            while(!out_queue.empty()){
                out_queue.pop(buff_idx);
            }

            in_queue.done();
            out_queue.done();

            close(fd);
            //std::cout << "File reading is ready\n";
        };

        //lambda functions that hash the kmers in a text chunk
        //note: it is not necessary for this function to be a lambda. It can be a static function

        // MODIFIED LONG
        auto bloom_filter_kmers =[&](chunk_type& chunk, size_t format){

            off_t i =0, last;
            size_t n_strings=0;

            // Rolling hasher for hash table positions
            //RollingHasherDual* rolling_hasher = new RollingHasherDual(ht_size, kmer_len);
            // Rolling hasher for bloom filter root hashes 
            //RollingHasherDual* bf_rolling_hasher = new RollingHasherDual(rolling_hasher_mod, kmer_len, bf_modmulinv, bf_multiplier, ht_size);
            RollingHasherDual* bf_rolling_hasher = new RollingHasherDual(rolling_hasher_mod, kmer_len, bf_modmulinv, bf_multiplier, ht_size, true);
            // Vector for storing hash values
            //std::vector<uint64_t> bloom_filter_hash_values(hash_functions, 0);
            std::vector<uint64_t> adbf_hash_values(hash_functions, 0);

            // --- Build k-mer factory ---
            KMerFactoryCanonical2BC* kmer_factory = new KMerFactoryCanonical2BC(k);

            switch (format) {
                case PLAIN://one-string-per-line format
                {
                    uint64_t new_char = 0;
                    assert(chunk.syms_in_buff>=k);
                    //slide a window over the buffer
                    while(i<chunk.syms_in_buff){
                        new_char =  uint64_t(twobitstringfunctions::char2int(chunk.buffer[i]));
                        
                        if (new_char > 3ULL){
                            kmer_factory->reset();
                        } else {
                            kmer_factory->push_new_integer(new_char);
                        }
                            
                        if (kmer_factory->get_number_of_stored_characters() == 0){
                            bf_rolling_hasher->reset();
                        } else {
                            bf_rolling_hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
                        }

                        if (kmer_factory->get_number_of_stored_characters() == kmer_len)
                        {
                            // Insert the k-mer in the bloom filter
                            uint64_t current_root_hash = std::min(bf_rolling_hasher->get_current_hash_backward_rqless(), bf_rolling_hasher->get_current_hash_forward_rqless());
                            bf->insertion_process(current_root_hash, adbf_hash_values);
                        }
                        i++;
                    }
                    if (print_other_stuff)
                        std::cout << "Chunk done\n";
                    break;
                }
                case FASTA: //fasta formta
                {
                    uint64_t new_char = 0;
                    assert(chunk.syms_in_buff>=k);
                    bool parsing_header = chunk.broken_header;    
                    //slide a window over the buffer
                    while(i<chunk.syms_in_buff){
                        // If the current character is header starting character, reset read buffer
                        if (chunk.buffer[i]=='>')
                        {
                            parsing_header = true;
                        }

                        if (parsing_header)
                        {
                            while((i<chunk.syms_in_buff) && (chunk.buffer[i]!='\n'))
                                i++;
                            i++;
                            parsing_header = false;
                            kmer_factory->reset();
                            bf_rolling_hasher->reset();
                            continue;
                        }                        
                        // If the next character is newline, skip it
                        if (chunk.buffer[i]=='\n')
                        {
                            i++;
                            continue;
                        }
                        new_char =  uint64_t(twobitstringfunctions::char2int(chunk.buffer[i]));
                        if (new_char > 3ULL){
                            kmer_factory->reset();
                        } else {
                            kmer_factory->push_new_integer(new_char);
                        }
                            
                        if (kmer_factory->get_number_of_stored_characters() == 0){
                            bf_rolling_hasher->reset();
                        } else {
                            bf_rolling_hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
                        }
                        
                        if (kmer_factory->get_number_of_stored_characters() == kmer_len)
                        {
                            // Insert the k-mer in the bloom filter
                            uint64_t current_root_hash = std::min(bf_rolling_hasher->get_current_hash_backward_rqless(), bf_rolling_hasher->get_current_hash_forward_rqless());
                            bf->insertion_process(current_root_hash, adbf_hash_values);
                        } 
                        
                        i++;
                    }
                    if (print_other_stuff)
                        std::cout << "Chunk " << chunk.id << " done\n";
                    break;
                }
                case FASTQ: //fastq format
                    //TODO
                    std::cout<<"Not implemented yet"<<std::endl;
                    break;
                default:
                    std::cout<<"Error : format not recognized"<<std::endl;
                    break;

            }
            delete kmer_factory;
            delete bf_rolling_hasher;
        };

        //lambda function that gets chunks from the IN queue and calls the hash_kmers lambda
        //we feed this function to std::thread
        auto string_worker = [&](size_t worker_id){

            size_t buff_id;
            bool res;
            size_t consumed_kmers = 0;

            while(true){
                res = in_queue.pop(buff_id);//the thread will wait until there is something to pop
                assert(text_chunks[buff_id].bytes>0);
                if(!res) break;
                bloom_filter_kmers(text_chunks[buff_id], format);
                consumed_kmers+=text_chunks[buff_id].syms_in_buff-k+1;
                out_queue.push(buff_id);//the thread will wait until the stack is free to push
            }

            /*
            if (print_other_stuff)
            {//TODO just testing
                std::unique_lock lck(mtx);
                std::cout<<"Thread "<<worker_id<<" consumed "<<consumed_kmers<<" kmers "<<std::endl;
            }
            */
        };

        std::vector<std::thread> threads;
        threads.emplace_back(io_worker);
        for(size_t i=0;i<n_threads;i++){
            threads.emplace_back(string_worker, i);
        }

        for(auto & thread : threads){
            //if (thread.joinable())
            thread.join();
        }

        std::vector<chunk_type>().swap(text_chunks);

        
        //remove the pages of the input file from the page cache
#ifdef __linux__
        posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
        close(fd);

        auto end_filtering = std::chrono::high_resolution_clock::now();

        if (print_times)
        {
            auto filtering_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_filtering - start_filtering);
            std::cout << "Time used to bloom filter k-mers: " << filtering_duration.count() << " microseconds\n";
        }    
    }
};

#endif //PARALLEL_PARSING_PARALLEL_PARSER_BF_HPP
