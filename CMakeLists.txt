cmake_minimum_required(VERSION 3.10)

# set project name
project(KaarmeKMerHashtable)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
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
        external/xxHash/xxhash.c
        )
        

target_link_libraries(kaarme pthread z)
target_compile_options(kaarme PRIVATE
        -Wall 
        -Wextra 
        -Wpedantic
        #-ggdb
        #-O1
        -O3
        #-ffast-math 
        #-march=native
        -Wno-unused-variable
        -Wno-unused-parameter
        #-flto
        #-funroll-loops
        #-finline-functions
        #-fno-stack-protector
        #-funswitch-loops
        #-fprefetch-loop-arrays
        )

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "Adding -march=native compiler flag due to ${CMAKE_CXX_COMPILER_ID}")
        target_compile_options(kaarme PUBLIC -march=native -Wno-unused-but-set-variable)

endif()

target_include_directories(kaarme PUBLIC include external/xxHash)
                           
