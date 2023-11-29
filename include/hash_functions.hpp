#include <vector>
#include <cstdint>
#include "kmer_factory.hpp"
#include "functions_math.hpp"
#include "functions_strings.hpp"

#pragma once

class RollingHasherDual2P
{
     public:

        // Hash table size (mod)
        uint64_t q;
        // Mod for return value
        uint64_t rq;
        // Current hash value, forward version
        __uint128_t current_hash_forward;
        // Current hash value, backward version        
        __uint128_t current_hash_backward;
        // k-mer length
        uint64_t m;
        // Bits per character
        uint64_t bpc;
        // Mask for one character
        uint64_t character_mask;
        // Alphabet size
        uint64_t d;
        // Multiplicative inverse of d under mod q
        __uint128_t di; 
        // Big Power
        __uint128_t h;
        // Hashed rollers
        uint64_t hashed_count;


        // Constructor
        RollingHasherDual2P(uint64_t q, uint64_t m);
        // Constructor 2
        RollingHasherDual2P(uint64_t q, uint64_t m, uint64_t modular_multiplicative_inverse, uint64_t multiplier);
        // Constructor 3
        RollingHasherDual2P(uint64_t q, uint64_t m, uint64_t modular_multiplicative_inverse, uint64_t multiplier, uint64_t return_q);
        // Destructor
        ~RollingHasherDual2P(){}
        // Update rolling hash with only the incoming character
        void update_rolling_hash_in(uint64_t in);
        // Update rolling hash with both incoming and outgoing characters
        void update_rolling_hash_in_and_out(uint64_t in, uint64_t out);
        // Update rolling hash (UNIVERSAL)
        void update_rolling_hash(uint64_t in, uint64_t out);
        // Return the current forward hash value
        uint64_t get_current_hash_forward();
        // Return the current backward hash value
        uint64_t get_current_hash_backward();
        // Return the current forward hash value
        uint64_t get_current_hash_forward_rqless();
        // Return the current backward hash value
        uint64_t get_current_hash_backward_rqless();
        // Reset hasher state
        void reset();
        // Load full contents from canonical k-mer factory
        void load_full_factory_canonical(KMerFactoryCanonical2BC* kmer_factory);
};

class RollingHasherDual
{
     public:

        // Hash table size (mod)
        uint64_t q;
        // Mod for return value
        uint64_t rq;
        // Current hash value, forward version
        __uint128_t current_hash_forward;
        // Current hash value, backward version        
        __uint128_t current_hash_backward;
        // k-mer length
        uint64_t m;
        // Bits per character
        uint64_t bpc;
        // Mask for one character
        uint64_t character_mask;
        // Alphabet size
        uint64_t d;
        // Multiplicative inverse of d under mod q
        __uint128_t di; 
        // Big Power
        __uint128_t h;
        // Hashed rollers
        uint64_t hashed_count;
        // Two bit mod
        bool tbm;


        // Constructor
        RollingHasherDual(uint64_t q, uint64_t m);
        // Constructor 2
        RollingHasherDual(uint64_t q, uint64_t m, uint64_t modular_multiplicative_inverse, uint64_t multiplier);
        // Constructor 3
        RollingHasherDual(uint64_t q, uint64_t m, uint64_t modular_multiplicative_inverse, uint64_t multiplier, uint64_t return_q);
        // Constructor 3
        RollingHasherDual(uint64_t q, uint64_t m, uint64_t modular_multiplicative_inverse, uint64_t multiplier, uint64_t return_q, bool twobitmod);
        // Destructor
        ~RollingHasherDual(){}
        // Update rolling hash with only the incoming character
        void update_rolling_hash_in(uint64_t in);
        // Update rolling hash with both incoming and outgoing characters
        void update_rolling_hash_in_and_out(uint64_t in, uint64_t out);
        // Update rolling hash (UNIVERSAL)
        void update_rolling_hash(uint64_t in, uint64_t out);
        // Return the current forward hash value
        uint64_t get_current_hash_forward();
        // Return the current backward hash value
        uint64_t get_current_hash_backward();
        // Return the current forward hash value
        uint64_t get_current_hash_forward_rqless();
        // Return the current backward hash value
        uint64_t get_current_hash_backward_rqless();
        // Reset hasher state
        void reset();
        // Load full contents from canonical k-mer factory
        void load_full_factory_canonical(KMerFactoryCanonical2BC* kmer_factory);
};

class AdderHasher1
{   
    private:

        uint64_t hash_table_slots;
    
    public:

        AdderHasher1(uint64_t slots);

        ~AdderHasher1(){};

        uint64_t calculate_hash(KMerFactoryCanonical2BC * kmer_factory);

};

class ProbeHasher1
{
    public:

        ProbeHasher1(){};

        ~ProbeHasher1(){};

        // Double hahsing probing
        uint64_t probe_1(uint64_t item, uint64_t M);

        // Quadratic probing
        uint64_t probe_2(uint64_t iteration);

        // Linear probing
        uint64_t probe_3(uint64_t iteration);

        // Better quadratic probing
        uint64_t probe_4(uint64_t iteration, uint64_t position, uint64_t modulo);

};