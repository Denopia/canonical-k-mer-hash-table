//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_PARALLEL_PARSER_HPP
#define PARALLEL_PARSING_PARALLEL_PARSER_HPP

#include "text_reader.h"
#include "ts_queue.h"
#include <thread>
#include <zlib.h>

enum bio_format{FASTA=0, FASTQ=1, PLAIN=2};

template<class sym_type,
         bool is_gzipped=false>
struct parse_input{

    void operator()(std::string& input_file, off_t chunk_size, size_t active_chunks, size_t n_threads, off_t k){

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
        //

        //get the file size
        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return;

        size_t format = PLAIN; //manage to get the input format
        std::mutex mtx; //just for debugging (you can remove it afterwards)

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

            while(chunk_id<active_chunks && rem_bytes>=k){

                tmp_ck_size = std::min(chunk_size, rem_bytes);
                text_chunks[chunk_id].bytes = tmp_ck_size;
                text_chunks[chunk_id].buffer = (sym_type *)malloc(tmp_ck_size);
                text_chunks[chunk_id].id = chunk_id;

                if constexpr (is_gzipped){
                    rem_bytes = read_chunk_from_gz_file<chunk_type>(gfd, text_chunks[chunk_id], rem_bytes, k-1);
                }else{
                    rem_bytes = read_chunk_from_file<chunk_type>(fd, text_chunks[chunk_id], rem_bytes, k-1);
                }
                in_queue.push(chunk_id);//as soon as we push, the chunks become visible of the worker threads to consume them
                chunk_id++;
            }

            size_t buff_idx;
            while(rem_bytes>=k){
                out_queue.pop(buff_idx);//it will wait until out_strings contains something
                text_chunks[buff_idx].id = chunk_id++;

                if constexpr (is_gzipped){
                    rem_bytes = read_chunk_from_gz_file<chunk_type>(gfd, text_chunks[buff_idx], rem_bytes, k-1);
                }else{
                    rem_bytes = read_chunk_from_file<chunk_type>(fd, text_chunks[buff_idx], rem_bytes, k-1);
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
        };

        //lambda functions that hash the kmers in a text chunk
        //note: it is not necessary for this function to be a lambda. It can be a static function
        auto hash_kmers =[&](chunk_type& chunk, size_t format){

            off_t i =0, last;

            switch (format) {
                case PLAIN://one-string-per-line format
                    assert(chunk.syms_in_buff>=k);

                    while(i<k){
                        //TODO compute the fingerprint
                        i++;
                    }

                    //slide a window over the buffer
                    while(i<chunk.syms_in_buff){

                        if(chunk.buffer[i]=='\n'){//we consumed a string
                            last = i+k; //rightmost position of the kmer in the prefix of the next string
                            for(off_t u=i+1;u<=last;u++){
                                 //TODO compute the fingerprint from scratch for the next string
                            }
                        }else{
                            //TODO update the fingerprint
                        }

                        //TODO hash the kmer
                        i++;
                    }
                    break;
                case FASTA: //fasta formta
                    //TODO
                    std::cout<<"Not implemented yet"<<std::endl;
                    break;
                case FASTQ: //fastq format
                    //TODO
                    std::cout<<"Not implemented yet"<<std::endl;
                    break;
                default:
                    std::cout<<"Error : format not recognized"<<std::endl;
                    break;

            }
        };

        //lambda function that gets chunks from the IN queue and calls the hash_kmers lambda
        //we feed this function to std::thread
        auto string_worker = [&](){

            size_t buff_id;
            bool res;
            //size_t consumed_kmers = 0;

            while(true){
                res = in_queue.pop(buff_id);//the thread will wait until there is something to pop
                assert(text_chunks[buff_id].bytes>0);
                if(!res) break;
                hash_kmers(text_chunks[buff_id], format);
                //consumed_kmers+=text_chunks[buff_id].syms_in_buff-k+1;
                out_queue.push(buff_id);//the thread will wait until the stack is free to push
            }

            /*{
                std::unique_lock lck(mtx);
                std::cout<<"This thread consumed "<<consumed_kmers<<" kmers "<<std::endl;
            }*/
        };

        std::vector<std::thread> threads;
        threads.emplace_back(io_worker);
        for(size_t i=0;i<n_threads;i++){
            threads.emplace_back(string_worker);
        }

        for(auto & thread : threads){
            thread.join();
        }

        std::vector<chunk_type>().swap(text_chunks);

        //remove the pages of the input file from the page cache
#ifdef __linux__
        posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif
        close(fd);
    }
};
#endif //PARALLEL_PARSING_PARALLEL_PARSER_HPP
