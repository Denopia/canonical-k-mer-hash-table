#include <vector>
#include <cstdint>
#include "kmer_factory.hpp"
#include "functions_math.hpp"
#include "functions_strings.hpp"

#pragma once


class RollingHasherDual
{
     public:

        // Hash table size (mod)
        uint64_t q;
        // Current hash value, forward version
        uint64_t current_hash_forward;
        // Current hash value, backward version        
        uint64_t current_hash_backward;
        // k-mer length
        uint64_t m;
        // Bits per character
        uint64_t bpc;
        // Mask for one character
        uint64_t character_mask;
        // Alphabet size
        uint64_t d;
        // Multiplicative inverse of d under mod q
        uint64_t di; 
        // Big Power
        uint64_t h;
        // Hashed rollers
        uint64_t hashed_count;

        // Constructor
        RollingHasherDual(uint64_t q, uint64_t m);
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