//
// Created by Diaz, Diego on 31.3.2023.
//

#ifndef PARALLEL_PARSING_TEXT_READER_H
#define PARALLEL_PARSING_TEXT_READER_H

#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <cassert>
#include <iostream>


template<class sym_t>
struct text_chunk {
    typedef sym_t sym_type;
    size_t id{}; //check of the ID (relative to the text)
    off_t bytes{}; //max. number of bytes that can be loaded in the chunk
    sym_t * buffer = nullptr; //buffer containing chunk data
    off_t syms_in_buff{}; //number of elements in the buffer
    bool broken_header=false;//bool indicating if in the previous chunk there was a line break after the rightmost header symbol

    ~text_chunk(){
        if(buffer!= nullptr){
            free(buffer);
        }
    }

    inline sym_type operator[](size_t idx) const {
        assert(idx<syms_in_buff);
        return buffer[idx];
    }
};

template<class text_chunk_t,
         typename sym_t = typename text_chunk_t::sym_type>
off_t read_chunk_from_gz_file(gzFile gfd, // file descriptor
                              text_chunk_t& chunk, // reference to the chunk struct were the data will be stored
                              off_t rem_text_bytes, // number of remaining bytes in the file
                              off_t k,
                              bool& next_has_broken_header,
                              sym_t start_symbol) { //reposition the file k positions back (for the kmers)

    off_t chunk_bytes = chunk.bytes<rem_text_bytes ? chunk.bytes : rem_text_bytes;
    chunk.bytes = chunk_bytes;

    off_t acc_bytes = 0;
    off_t read_bytes;
    off_t fd_buff_bytes = 8388608;// 8MiB buffer

    off_t sym_bytes = sizeof(sym_t);
    off_t buff_pos=0;

    while(chunk_bytes>0) {

        fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
        read_bytes = gzread(gfd, &chunk.buffer[buff_pos], fd_buff_bytes);
        assert(read_bytes>0);

        buff_pos = buff_pos + (read_bytes/sym_bytes);
        chunk_bytes-=read_bytes;
        acc_bytes+=read_bytes;
    }

    chunk.syms_in_buff =  acc_bytes/sym_bytes;

    if(start_symbol!=0){
        next_has_broken_header = true;
        size_t i = chunk.syms_in_buff;
        while(i-->0 && chunk.buffer[i]!=start_symbol){
            if(chunk.buffer[i]=='\n'){
                next_has_broken_header = false;
                break;
            }
        }
    }

    gzseek(gfd, k*-1, SEEK_CUR);
    acc_bytes-=k;
    rem_text_bytes-=acc_bytes;

    //std::cout<<acc_bytes<<" "<<rem_text_bytes<<std::endl;
    assert(chunk_bytes==0);

    return rem_text_bytes;
}

template<class text_chunk_t,
         typename sym_t = typename text_chunk_t::sym_type>
off_t read_chunk_from_file(int fd, // file descriptor
                          text_chunk_t& chunk, // reference to the chunk struct were the data will be stored
                          off_t rem_text_bytes, // number of remaining bytes in the file
                          off_t k, //reposition the file description k-1 positions back (for the kmers)
                          bool &next_has_broken_header,
                          sym_t start_symbol=0) {

    //std::cout << "Start preparing text chunk\n";
    //std::cout << "Remaining bytes = " << rem_text_bytes << "\n";
    off_t chunk_bytes = chunk.bytes<rem_text_bytes ? chunk.bytes : rem_text_bytes;
    chunk.bytes = chunk_bytes;

    off_t acc_bytes = 0;
    off_t read_bytes;
    off_t fd_buff_bytes = 8388608;// 8MiB buffer

    off_t sym_bytes = sizeof(sym_t);
    off_t buff_pos=0;

    while(chunk_bytes>0) {

        fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
        read_bytes = read(fd, &chunk.buffer[buff_pos], fd_buff_bytes);
        assert(read_bytes>0);

        buff_pos = buff_pos + (read_bytes/sym_bytes);
        chunk_bytes-=read_bytes;
        acc_bytes+=read_bytes;
    }
    
    //std::cout << "----------------------------------------------------------------\n";
    //std::cout << "acc bytes 1 = " << acc_bytes << "\n";
    //std::cout << "chunk bytes bytes = " << chunk_bytes << "\n";

    chunk.syms_in_buff =  acc_bytes/sym_bytes;

    off_t fake_symbols = 0;
    off_t real_symbols = 0;
    off_t si = chunk.syms_in_buff-1;

    //std::cout << "FAKE SYMBOLS = " << fake_symbols << "\n";

    if(start_symbol!=0){
        // First, check how many "fake symbols" i.e. newline symbols are encountered 
        // starting from the end of the text chunk until k-1 "real" symbols are seen.
        // This is only relevant for fasta files since they allow newline characters in the middle of a read sequence(?)
        if (start_symbol == '>')
        {
            while(real_symbols < k && si >= 0)
            {
                if (chunk.buffer[si] != '\n')
                    real_symbols++;
                else
                    fake_symbols++;
                si--;
            }
        }

        // Chunk does not have enough real symbols
        if (real_symbols != k)
        {
            rem_text_bytes = 0;
            return rem_text_bytes;
        }
        
        // Next, check if the text chunk ends with a header.
        // Skip the last k-1 characters to avoid weird behavior(??)
        next_has_broken_header = true;
        size_t i = chunk.syms_in_buff - 1 - k - fake_symbols; 
        //size_t i = chunk.syms_in_buff;
        //while(i-->0 && chunk.buffer[i]!=start_symbol){
        while(i >= 0) {
            if(chunk.buffer[i]==start_symbol)
                break;
            if(chunk.buffer[i]=='\n'){
                //std::cout << "Newline found, not broken header\n";
                next_has_broken_header = false;
                break;
            }
            i--;
        }
    } else {
        //std::cout << "No header checking\n";
    }

    // If the number of characters in the text chunk minus newline characters 
    // is less than k the chunk does not contain a full k-mer
    //if (acc_bytes-fake_symbols < k+1){
    //    rem_text_bytes = 0;
    //    return rem_text_bytes;
    //}

    //if (next_has_broken_header)
    //    std::cout << "\nHEADER IS BROKEN\n\n";

    //lseek(fd, k*-1, SEEK_CUR);
    // seek backward k-1 characters and the number of fake symbols
    lseek(fd, ((-1*k)-fake_symbols), SEEK_CUR);
    //acc_bytes-=k;
    // acc_bytes(?) is also decresed by the number of fake symbols
    acc_bytes= acc_bytes - k - fake_symbols;
    //std::cout << "acc bytes 2 = " << acc_bytes << "\n";
    //std::cout << "rem text bytes = " << rem_text_bytes << "\n";
    
    if (acc_bytes == 0)
        rem_text_bytes = 0;
    else
        rem_text_bytes-=acc_bytes;

    assert(chunk_bytes==0);

    //std::cout << "Chunk prepared successfully\n";
    return rem_text_bytes;
}
#endif //PARALLEL_PARSING_TEXT_READER_H
