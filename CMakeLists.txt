cmake_minimum_required(VERSION 3.10)

# set project name
project(KMerHashtable)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(EXTRA_INCLUDES include)

add_executable(kht main.cpp
        source/file_reader.cpp
        source/hash_functions.cpp
        source/functions_strings.cpp
        source/kmer_factory.cpp
        source/functions_math.cpp
        source/program_runs.cpp
        source/kmer.cpp
        source/kmer_hash_table.cpp)

target_link_libraries(kht pthread z)
target_compile_options(kht PRIVATE
        -Wall -Wextra -Wpedantic
        -O3
        -ffast-math -ggdb)
target_include_directories(kht PUBLIC include)
                           
