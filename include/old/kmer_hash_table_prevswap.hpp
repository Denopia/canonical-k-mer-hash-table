#include <vector>
#include <cstdint>
#include "kmer_factory.hpp"
#include "vectors.hpp"
#include "kmer.hpp"
#include "hash_functions.hpp"
#include "functions_math.hpp"
#include "bit_vectors.hpp"

#pragma once

class BasicHashTable1
{

    private:
        // Size of the hash table
        uint32_t size;
        // Length of the k-mers
        uint64_t kmer_len;
        // Blocks taken by one k-mer
        uint64_t single_kmer_blocks;
        // Array that holds the hash table k-mers
        uint64_t* hash_table_array_kmers;
        //sdsl::int_vector<64> hash_table_array_kmers;
        // Array that holds the hash table k-mer counts
        uint32_t* hash_table_array_counts;
        // Number of bits used to represent one character
        uint64_t bits_per_char;
        // Number of items inserted in the hash table
        uint64_t inserted_items;
        // Number of complete k-mers inserted in the hash table (subset of inserted items)
        uint64_t inserted_complete_kmers;
        // k-mers with more than 1 occurrence
        uint64_t solid_kmers;

    public:

        BasicHashTable1(uint64_t s, uint64_t k);

        ~BasicHashTable1();

        uint32_t get_kmer_count_in_slot(uint64_t slot);

        void resize();

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_complete_kmers();

        uint64_t get_number_of_solid_kmers();

        uint64_t insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);

        uint64_t find(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);

        uint64_t find_and_increment(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);

        uint64_t find_and_increment_andifnot_insert(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);
};

class PointerHashTable1
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMer* hash_table_array;
        // Size of the hash table
        uint32_t size;
        // Length of the k-mers
        uint64_t kmer_len;
        // Number of bits used to represent one character
        uint64_t bits_per_char;
        // Number of items inserted in the hash table
        uint64_t inserted_items;
        // Number of complete k-mers inserted in the hash table (subset of inserted items)
        uint64_t inserted_complete_kmers;
        // Number of solid k-mer (count at least 2)
        uint64_t solid_kmers;


        ProbeHasher1 * probe_hasher;
        uint64_t probing_prime;

    public:

        PointerHashTable1(uint64_t s, uint64_t k);

        ~PointerHashTable1();

        uint32_t get_kmer_count_in_slot(uint64_t slot);

        void resize();

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_complete_kmers();

        uint64_t get_number_of_solid_kmers();

        // Check if the currently loaded k-mer is in the given slot
        bool kmer_slot_check(KMerFactory2BC* kmer_factory, uint64_t slot);

        // Find where a k-mer is in the hash table
        uint64_t find(KMerFactory2BC* kmer_factory, uint64_t kmer_hash);

        // Find where a k-mer is in the hash table
        uint64_t find_using_prev(uint32_t new_chracter, uint64_t kmer_slot, uint64_t previous_kmer_slot);

        // Insert new k-mer to the hash table
        uint64_t insert_new_kmer(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);

        // Find the slot of the given k-mer, and increment count if found
        // and if not found -> insert as new k-mer
        uint64_t find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);
};

class PointerHashTable2
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMer* hash_table_array;
        // Size of the hash table
        uint32_t size;
        // Length of the k-mers
        uint64_t kmer_len;
        // Number of bits used to represent one character
        uint64_t bits_per_char;
        // Number of items inserted in the hash table
        uint64_t inserted_items;
        // Number of complete k-mers inserted in the hash table (subset of inserted items)
        uint64_t inserted_complete_kmers;
        // Number of solid k-mer (count at least 2)
        uint64_t solid_kmers;

        // Probing related stuff
        ProbeHasher1 * probe_hasher;
        uint64_t probing_prime;

        // Temp array related stuff
        uint32_t max_temp_slots;
        uint32_t touched_temp_slots;
        uint32_t temp_slots_in_use;
        uint32_t max_temp_slot_in_use;
        uint32_t smallest_unused_temp_slot;
        std::vector<uint64_t> temp_array;
        sdsl::bit_vector temp_free_slots;

    public:

        PointerHashTable2(uint64_t s, uint64_t k, uint64_t b);

        ~PointerHashTable2();

        uint32_t get_kmer_count_in_slot(uint64_t slot);

        void resize();

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_complete_kmers();

        uint64_t get_number_of_solid_kmers();

        uint64_t get_number_of_max_temp_slots();

        uint64_t get_number_of_touched_temp_slots();

        uint64_t get_number_of_temp_slots_in_use();

        uint64_t get_number_of_max_temp_slots_in_use();

        uint64_t get_smallest_unused_temp_slot();

        // Check if the currently loaded k-mer is in the given slot
        int kmer_slot_check(KMerFactory2BC* kmer_factory, uint64_t slot);

        // Find where a k-mer is in the hash table
        uint64_t find(KMerFactory2BC* kmer_factory, uint64_t kmer_hash);
        
        // Find where a k-mer is in the hash table
        uint64_t find_using_prev(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);

        // Find where a k-mer is in the hash table
        uint64_t find_and_modify(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);

        // Insert new k-mer to the hash table main array
        uint64_t insert_new_kmer_in_main(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);
        
        // Insert new item to the temp array
        uint64_t insert_new_kmer_in_temp(KMerFactory2BC* kmer_factory, uint64_t kmer_hash);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);

        // Find the slot of the given k-mer, and increment count if found
        // and if not found -> insert as new k-mer
        uint64_t find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);

        void clean_up_temp_array(KMerFactory2BC* kmer_factory, RollingHasher1* hasher);

        char get_oldest_char_in_temp_slot_and_shift_left(uint32_t slot, KMerFactory2BC* kmer_factory);
};
