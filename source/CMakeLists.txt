cmake_minimum_required(VERSION 3.10)

# set project name
project(KMerHashtable)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include_directories(../sdsl-lite-lib/include/sdsl ../sdsl-lite-lib/include)
#find_library(SDSL )
#set(EXTRA_LIBS ${PROJECCT_SOURCE_DIR}/../sdsl-lite-lib/lib/libsdsl)
#add_library(sdsl divsufsort divsufsort64 )
#file(GLOB SOURCES_FILES "source/*.cpp" main.cpp)

#message(STATUS "apua1")
#message(STATUS "apua2")

# add executable

set(EXTRA_INCLUDES ../include)

add_executable(${PROJECT_NAME} main.cpp file_reader.cpp hash_functions.cpp functions_strings.cpp kmer_factory.cpp functions_math.cpp program_runs.cpp kmer.cpp kmer_hash_table.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC ${EXTRA_LIBS} ${LIBSDSL_LIBRARIES})

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -O3 -ffast-math -ggdb)

target_include_directories(KMerHashtable PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           ${EXTRA_INCLUDES}
                           ${LIBSDSL_INCLUDE_DIRS}
                           )
                           
