cmake_minimum_required(VERSION 3.10)

# set project name
project(KaarmeKMerHashtable)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(EXTRA_INCLUDES include)

add_executable(kaarme main.cpp
        source/file_reader.cpp
        source/hash_functions.cpp
        source/functions_strings.cpp
        source/kmer_factory.cpp
        source/functions_math.cpp
        source/functions_kmer_mod.cpp
        source/functions_bloom_filter.cpp
        source/kmer.cpp
        source/kmer_hash_table.cpp
        external/xxhash.c
        )
        

target_link_libraries(kaarme pthread z)
target_compile_options(kaarme PRIVATE
        -Wall 
        -Wextra 
        -Wpedantic
        -ggdb
        #-O1
        -O3
        #-ffast-math 
        -march=native
        #-flto
        #-funroll-loops
        #-finline-functions
        #-fno-stack-protector
        #-funswitch-loops
        #-fprefetch-loop-arrays
        )
target_include_directories(kaarme PUBLIC include xxHash)
                           
