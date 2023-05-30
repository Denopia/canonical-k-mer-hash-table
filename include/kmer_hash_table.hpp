#include <vector>
#include <cstdint>
#include <fstream>
#include "kmer_factory.hpp"
//#include "vectors.hpp"
#include "kmer.hpp"
#include "hash_functions.hpp"
#include "functions_math.hpp"
//#include "bit_vectors.hpp"
#include <tuple>
//#include <sdsl/bit_vectors.hpp>
#include <atomic>
#include <bitset>

#pragma once


class BasicAtomicHashTable
{
    public:
        BasicAtomicKMer* kmer_array;
        uint64_t size;
        uint64_t kmer_len;
    
        BasicAtomicHashTable(uint64_t s, uint64_t k);

        ~BasicAtomicHashTable();

        //void insert_new_atomically(uint64_t kmer);

        //void update_count_atomically(uint64_t slot, uint64_t new_count);

        //uint64_t get_count(uint64_t slot);

        //uint64_t get_kmer(uint64_t slot);
};

class BasicAtomicFlagHashTableLong
{
    public:
        uint8_t* kmer_array;
        uint16_t* counts;
        std::atomic_flag* kmer_locks;
        uint64_t size;
        uint32_t kmer_len;
        uint32_t kmer_bytes;

        // s = size, k = k-mer length
        BasicAtomicFlagHashTableLong(uint64_t s, uint32_t k);

        ~BasicAtomicFlagHashTableLong();

        void write_kmers(uint64_t min_abundance, std::string& output_path);

        //void insert_new_atomically(uint64_t kmer);

        //void update_count_atomically(uint64_t slot, uint64_t new_count);

        //uint64_t get_count(uint64_t slot);

        //uint64_t get_kmer(uint64_t slot);
};

class BasicAtomicVariableHashTableLong
{
    public:
        uint8_t* kmer_array;
        std::atomic<uint32_t> * counts;
        uint64_t size;
        uint32_t kmer_len;
        uint32_t kmer_bytes;

        // s = size, k = k-mer length
        BasicAtomicVariableHashTableLong(uint64_t s, uint32_t k);

        ~BasicAtomicVariableHashTableLong();

        void write_kmers(uint64_t min_abundance, std::string& output_path);

        //void insert_new_atomically(uint64_t kmer);

        //void update_count_atomically(uint64_t slot, uint64_t new_count);

        //uint64_t get_count(uint64_t slot);

        //uint64_t get_kmer(uint64_t slot);
};


class BasicAtomicVariableHashTableLong64
{
    public:
        uint64_t* kmer_array;
        std::atomic<uint32_t> * counts;
        uint64_t size;
        uint32_t kmer_len;
        uint32_t kmer_blocks;

        // s = size, k = k-mer length
        BasicAtomicVariableHashTableLong64(uint64_t s, uint32_t k);

        ~BasicAtomicVariableHashTableLong64();

        //void insert_new_atomically(uint64_t kmer);

        //void update_count_atomically(uint64_t slot, uint64_t new_count);

        //uint64_t get_count(uint64_t slot);

        //uint64_t get_kmer(uint64_t slot);
};


class PointerHashTableCanonicalAF
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMerAtomicFlag* hash_table_array;
        // Size of the hash table
        uint64_t size;
        // Length of the k-mers
        uint64_t kmer_len;
        // Number of bits used to represent one character
        uint64_t bits_per_char;
        // Number of items inserted in the hash table
        uint64_t inserted_items;
        // Number of complete k-mers inserted in the hash table (subset of inserted items)
        uint64_t inserted_complete_kmers;
        // Integers needed to store full k-mer in 2bits per char representation
        uint64_t kmer_blocks;
        // Probing related stuff
        ProbeHasher1 * probe_hasher;
        uint64_t probing_prime;
        // Secondary array stuff
        uint64_t max_secondary_slots;
        uint64_t touched_secondary_slots;
        uint64_t secondary_slots_in_use;
        uint64_t max_secondary_slot_in_use;
        uint64_t smallest_unused_secondary_slot;
        std::vector<uint64_t> secondary_array;
        std::vector<uint8_t> secondary_free_slots;

        uint64_t max_kmer_reconstruction_chain;
        uint64_t total_reconstruction_chain;

        //std::atomic_flag* main_locks;
        //std::vector<std::atomic_flag> secondary_locks;
        std::atomic_flag secondary_lock;

    public:
        // s = slots, k = k-mer length, b = 64bit blocks per k-mer
        PointerHashTableCanonicalAF(uint64_t s, uint64_t k, uint64_t b);

        ~PointerHashTableCanonicalAF();

        uint64_t get_kmer_count_in_slot(uint64_t slot);

        bool kmer_in_slot_is_complete(uint64_t slot);

        bool slot_is_occupied(uint64_t slot);

        void resize();

        //uint64_t get_size();

        uint64_t process_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        // Find where a k-mer is in the hash table
        uint64_t find(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        int quick_kmer_slot_check_sus(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot, uint64_t predecessor_slot);

        bool full_kmer_slot_check_NO_SECONDARY(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot);

        bool full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot);

        uint64_t insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        void write_kmers_on_disk_separately(uint64_t min_abundance, std::string& output_path);

        void write_kmers_on_disk_in_blocks(KMerFactoryCanonical2BC* kmer_factory, uint64_t min_abundance, std::string& output_path);

        uint64_t insert_new_kmer_in_secondary(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher);

        uint64_t get_secondary_array_char(uint64_t secondary_array_position, int char_position);

        std::string reconstruct_kmer_in_slot(uint64_t slot);

        bool check_for_cycle(uint64_t reconstruction_slot, uint64_t avoid_slot);

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_items_in_main();

        uint64_t get_number_of_max_secondary_slots();

        uint64_t get_number_of_secondary_slots_in_use();

        uint64_t get_max_number_of_secondary_slots_in_use();

        void write_kmers_on_disk_separately_faster(uint64_t min_abundance, std::string& output_path);

        void write_kmers_on_disk_separately_even_faster(uint64_t min_abundance, std::string& output_path);

};


class PointerHashTableCanonicalAV
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMerAtomicVariable* hash_table_array;
        // Size of the hash table
        uint64_t size;
        // Length of the k-mers
        uint64_t kmer_len;
        // Number of bits used to represent one character
        uint64_t bits_per_char;
        // Number of items inserted in the hash table
        uint64_t inserted_items;
        // Number of complete k-mers inserted in the hash table (subset of inserted items)
        uint64_t inserted_complete_kmers;
        // Integers needed to store full k-mer in 2bits per char representation
        uint64_t kmer_blocks;
        // Probing related stuff
        ProbeHasher1 * probe_hasher;
        uint64_t probing_prime;
        // Secondary array stuff
        uint64_t max_secondary_slots;
        uint64_t touched_secondary_slots;
        uint64_t secondary_slots_in_use;
        uint64_t max_secondary_slot_in_use;
        uint64_t smallest_unused_secondary_slot;
        std::vector<uint64_t> secondary_array;
        std::vector<uint8_t> secondary_free_slots;

        uint64_t max_kmer_reconstruction_chain;
        uint64_t total_reconstruction_chain;

        std::atomic_flag secondary_lock;


    public:
        // s = slots, k = k-mer length, b = 64bit blocks per k-mer
        PointerHashTableCanonicalAV(uint64_t s, uint64_t k, uint64_t b);

        ~PointerHashTableCanonicalAV();

        uint64_t get_kmer_count_in_slot(uint64_t slot);

        bool kmer_in_slot_is_complete(uint64_t slot);

        bool slot_is_occupied(uint64_t slot);

        void resize();

        //uint64_t get_size();

        uint64_t process_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        // Find where a k-mer is in the hash table
        uint64_t find(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        int quick_kmer_slot_check_sus(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot, uint64_t predecessor_slot);

        bool full_kmer_slot_check_NO_SECONDARY(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot);

        bool full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot);

        uint64_t insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        void write_kmers_on_disk_separately(uint64_t min_abundance, std::string& output_path);

        void write_kmers_on_disk_in_blocks(KMerFactoryCanonical2BC* kmer_factory, uint64_t min_abundance, std::string& output_path);

        uint64_t insert_new_kmer_in_secondary(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher);

        uint64_t get_secondary_array_char(uint64_t secondary_array_position, int char_position);

        std::string reconstruct_kmer_in_slot(uint64_t slot);

        bool check_for_cycle(uint64_t reconstruction_slot, uint64_t avoid_slot);

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_items_in_main();

        uint64_t get_number_of_max_secondary_slots();

        uint64_t get_number_of_secondary_slots_in_use();

        uint64_t get_max_number_of_secondary_slots_in_use();

        void write_kmers_on_disk_separately_faster(uint64_t min_abundance, std::string& output_path);

        void write_kmers_on_disk_separately_even_faster(uint64_t min_abundance, std::string& output_path);

        uint64_t process_kmer_MT(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot);

        void analyze_pointer_chain_lengths();

        uint64_t count_reconstruction_chain_length_in_slot(uint64_t slot);

};



