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

template<class sym_t>
struct text_chunk {
    typedef sym_t sym_type;
    size_t id{}; //check of the ID (relative to the text)
    off_t bytes{}; //max. number of bytes that can be loaded in the chunk
    sym_t * buffer = nullptr; //buffer containing chunk data
    off_t syms_in_buff{}; //number of elements in the buffer

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

template<class text_chunk_t>
off_t read_chunk_from_gz_file(gzFile gfd, // file descriptor
                              text_chunk_t& chunk, // reference to the chunk struct were the data will be stored
                              off_t rem_text_bytes, // number of remaining bytes in the file
                              off_t k) { //reposition the file k positions back (for the kmers)

    using sym_t = typename text_chunk_t::sym_type;

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

    gzseek(gfd, k*-1, SEEK_CUR);
    acc_bytes-=k;
    rem_text_bytes-=acc_bytes;

    //std::cout<<acc_bytes<<" "<<rem_text_bytes<<std::endl;
    assert(chunk_bytes==0);

    return rem_text_bytes;
}

template<class text_chunk_t>
off_t read_chunk_from_file(int fd, // file descriptor
                          text_chunk_t& chunk, // reference to the chunk struct were the data will be stored
                          off_t rem_text_bytes, // number of remaining bytes in the file
                          off_t k) { //reposition the file description k-1 positions back (for the kmers)

    using sym_t = typename text_chunk_t::sym_type;

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

    chunk.syms_in_buff =  acc_bytes/sym_bytes;
    lseek(fd, k*-1, SEEK_CUR);
    acc_bytes-=k;
    rem_text_bytes-=acc_bytes;

    //std::cout<<acc_bytes<<" "<<rem_text_bytes<<" "<<chunk.syms_in_buff<<" "<<sym_bytes<<std::endl;

    assert(chunk_bytes==0);

    return rem_text_bytes;
}
#endif //PARALLEL_PARSING_TEXT_READER_H
