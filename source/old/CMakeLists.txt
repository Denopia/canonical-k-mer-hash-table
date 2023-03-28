cmake_minimum_required(VERSION 3.10)

# set project name
project(KMerHashtable)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

##################################################################################################################################
# SDSL-lite finder starts
##################################################################################################################################
# - Try to find LibSDSL
# Once done this will define
#  LIBSDSL_FOUND - System has SDSL
#  LIBSDSL_INCLUDE_DIRS - The SDSL include directories
#  LIBSDSL_LIBRARIES - The libraries needed to use SDSL
#  LIBSDSL_DEFINITIONS - Compiler switches required for using SDSL

find_package(PkgConfig)

pkg_check_modules(PC_LIBSDSL QUIET libsdsl)

set(LIBSDSL_DEFINITIONS ${PC_LIBSDSL_CFLAGS_OTHER})

find_path(LIBSDSL_INCLUDE_DIR sdsl/csa_sada.hpp
          HINTS ${PC_LIBSDSL_INCLUDEDIR}
                ${PC_LIBSDSL_INCLUDE_DIRS}
          PATHS
          		${PROJECT_SOURCE_DIR}/../sdsl-lite/include/
                #$ENV{HOME}/usr/include
                #$ENV{HOME}/include
          PATH_SUFFIXES libsdsl
          )

find_library(LIBSDSL_LIBRARY NAMES sdsl
		HINTS 	${PC_LIBSDSL_LIBDIR}
				${PC_LIBSDSL_LIBRARY_DIRS}
		PATHS
             	${PROJECT_SOURCE_DIR}/../sdsl-lite/lib
                #$ENV{HOME}/usr/lib
                #$ENV{HOME}/lib
                )

find_library(LIBDIVSUF_LIBRARY NAMES divsufsort
        HINTS ${PC_LIBSDSL_LIBDIR}
              ${PC_LIBSDSL_LIBRARY_DIRS}
        PATHS
        	  ${PROJECT_SOURCE_DIR}/../sdsl-lite/lib
              #$ENV{HOME}/usr/lib
              #$ENV{HOME}/lib
              )

find_library(LIBDIVSUF64_LIBRARY NAMES divsufsort64
        HINTS ${PC_LIBSDSL_LIBDIR}
        	  ${PC_LIBSDSL_LIBRARY_DIRS}
        PATHS
        	  ${PROJECT_SOURCE_DIR}/../sdsl-lite/lib
        	  #$ENV{HOME}/usr/lib
        	  #$ENV{HOME}/lib
        	  )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBSDSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibSDSL DEFAULT_MSG
        LIBSDSL_LIBRARY LIBDIVSUF_LIBRARY LIBDIVSUF64_LIBRARY LIBSDSL_INCLUDE_DIR)

mark_as_advanced(LIBSDSL_LIBRARY LIBDIVSUF_LIBRARY LIBDIVSUF64_LIBRARY LIBSDSL_INCLUDE_DIR)

set(LIBSDSL_LIBRARIES ${LIBSDSL_LIBRARY} ${LIBDIVSUF_LIBRARY} ${LIBDIVSUF64_LIBRARY})
set(LIBSDSL_INCLUDE_DIRS ${LIBSDSL_INCLUDE_DIR})

##################################################################################################################################
# SDSL-lite finder ends
##################################################################################################################################

#include_directories(../sdsl-lite-lib/include/sdsl ../sdsl-lite-lib/include)
#find_library(SDSL )
#set(EXTRA_LIBS ${PROJECCT_SOURCE_DIR}/../sdsl-lite-lib/lib/libsdsl)
#add_library(sdsl divsufsort divsufsort64 )
#file(GLOB SOURCES_FILES "source/*.cpp" main.cpp)

#message(STATUS "apua1")
message(STATUS "${LIBSDSL_INCLUDE_DIRS}")
#message(STATUS "apua2")

set(EXTRA_INCLUDES ../include ${PROJECT_BINARY_DIR}/../sdsl-lite/include/sdsl)

# add executable
add_executable(${PROJECT_NAME} main.cpp file_reader.cpp hash_functions.cpp kmer_factory.cpp kmer_hash_table.cpp functions_strings.cpp functions_math.cpp program_runs.cpp kmer.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC ${EXTRA_LIBS} ${LIBSDSL_LIBRARIES})

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -O3 -ffast-math)

target_include_directories(KMerHashtable PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           ${EXTRA_INCLUDES}
                           ${LIBSDSL_INCLUDE_DIRS}
                           )
                           
