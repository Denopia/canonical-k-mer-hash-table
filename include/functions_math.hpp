#include <cmath>
#include <cstdint>
#include <iostream>

#pragma once

//////////////////////////////////////////////
//
// This file has some useful math functions
//
//////////////////////////////////////////////

namespace mathfunctions
{
    /*
        This function returns the smallest prime 
        that is at least as large as the given integer

        Useful for finding hash table size 
        (maybe a prime at elast 1.33 times the expected amount of stored items)
    */
    uint64_t next_prime(uint64_t at_least);

    /*
        This function returns the smallest prime 
        that is at least as large as the given integer
        AND
        prime % 4 = 3

        Useful for finding hash table size 
        (maybe a prime at elast 1.33 times the expected amount of stored items)
    */
    uint64_t next_prime3mod4(uint64_t at_least);


    /*
        This function returns the multiplicative inverse of A
        under modulo M. 

        NOTE: M must be prime.
    */
    uint64_t modular_multiplicative_inverse(uint64_t A, uint64_t M);


    /*
        This function return x^y mod M
    */
    uint64_t power_under_modulo(uint64_t x, uint64_t y, uint64_t M);

    /*
        This function return the greatest common divisor between A and B
    */
    uint64_t greatest_common_divisor(uint64_t A, uint64_t B);

}