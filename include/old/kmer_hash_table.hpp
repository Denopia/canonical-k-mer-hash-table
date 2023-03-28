#include <vector>
#include <cstdint>
#include "kmer_factory.hpp"
#include "vectors.hpp"
#include "kmer.hpp"
#include "hash_functions.hpp"
#include "functions_math.hpp"
#include "bit_vectors.hpp"
#include <tuple>
#include <sdsl/bit_vectors.hpp>


#pragma once

class PointerHashTable2
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMer* hash_table_array;
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
        // Number of solid k-mer (count at least 2)
        uint64_t solid_kmers;

        uint64_t kmer_blocks;

        // Probing related stuff
        ProbeHasher1 * probe_hasher;
        uint64_t probing_prime;

        // Temp array related stuff
        uint64_t max_temp_slots;
        uint64_t touched_temp_slots;
        uint64_t temp_slots_in_use;
        uint64_t max_temp_slot_in_use;
        uint64_t smallest_unused_temp_slot;
        std::vector<uint64_t> temp_array;
        std::vector<uint8_t> temp_free_slots;
        //sdsl::bit_vector temp_free_slots;

    public:

        PointerHashTable2(uint64_t s, uint64_t k, uint64_t b);

        ~PointerHashTable2();

        uint64_t get_kmer_count_in_slot(uint64_t slot);

        uint64_t get_size();

        bool kmer_in_slot_is_complete(uint64_t slot);

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
        
        //Find where a k-mer is in the hash table using previous pointer to speedup
        std::tuple<bool,uint64_t> find_and_modify_using_prev(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);

        // Insert new k-mer to the hash table main array
        uint64_t insert_new_kmer_in_main(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);
        
        // Insert new item to the temp array
        uint64_t insert_new_kmer_in_temp(KMerFactory2BC* kmer_factory, uint64_t kmer_hash);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_ip_step);

        // Find the slot of the given k-mer, and increment count if found
        // and if not found -> insert as new k-mer
        uint64_t find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_up_step);

        void clean_up_temp_array(KMerFactory2BC* kmer_factory, RollingHasher1* hasher);

        void clean_up_temp_array_OLD(KMerFactory2BC* kmer_factory, RollingHasher1* hasher);

        char get_oldest_char_in_temp_slot_and_shift_left(uint64_t slot, KMerFactory2BC* kmer_factory);

        void write_kmers_on_disk(KMerFactory2BC* kmer_factory, uint64_t min_abundance, bool only_canonical, std::string output_path);

        //void copy_content_from_another_hash_table_ALSO_BROKEN(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table);
        
        //void copy_content_from_another_hash_table_BROKEN(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table);

        void copy_content_from_another_hash_table(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table);

        void copy_content_from_another_hash_table_OLD_AND_WORKING(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table);

        void copy_content_from_another_hash_table_BIGMEM(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table);

};


class PointerHashTableCanonical
{

    private:

        // k-mers are stored in this array
        OneCharacterAndPointerKMer* hash_table_array;
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
        // Temp array related stuff
        uint64_t max_temp_slots;
        uint64_t touched_temp_slots;
        uint64_t temp_slots_in_use;
        uint64_t max_temp_slot_in_use;
        uint64_t smallest_unused_temp_slot;
        std::vector<uint64_t> temp_array;
        std::vector<uint8_t> temp_free_slots;

    public:

        PointerHashTableCanonical(uint64_t s, uint64_t k, uint64_t b);

        ~PointerHashTableCanonical();

        uint64_t get_kmer_count_in_slot(uint64_t slot);

        uint64_t get_size();

        bool kmer_in_slot_is_complete(uint64_t slot);

        void resize();

        uint64_t get_number_of_inserted_items();

        uint64_t get_number_of_inserted_complete_kmers();

        uint64_t get_number_of_max_temp_slots();

        uint64_t get_number_of_touched_temp_slots();

        uint64_t get_number_of_temp_slots_in_use();

        uint64_t get_number_of_max_temp_slots_in_use();

        uint64_t get_smallest_unused_temp_slot();

        // Check if the currently loaded k-mer is in the given slot
        int full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t slot);

        // Find where a k-mer is in the hash table
        uint64_t find(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);
        
        // Find where a k-mer is in the hash table
        uint64_t find_using_prev(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);

        // Find where a k-mer is in the hash table
        uint64_t find_and_modify(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);
        
        //Find where a k-mer is in the hash table using previous pointer to speedup
        std::tuple<bool,uint64_t> find_and_modify_using_prev(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot);

        // Insert new k-mer to the hash table main array
        uint64_t insert_new_kmer_in_main(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot);
        
        // Insert new item to the temp array
        uint64_t insert_new_kmer_in_temp(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash);

        // Find the slot of the given k-mer, and increment count if found
        uint64_t find_and_increment(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_ip_step);

        // Find the slot of the given k-mer, and increment count if found
        // and if not found -> insert as new k-mer
        uint64_t add_kmer(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_up_step);

        void clean_up_temp_array(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher);

        char get_oldest_char_in_temp_slot_and_shift_left(uint64_t slot, KMerFactoryCanonical2BC* kmer_factory);

        void write_kmers_on_disk(KMerFactoryCanonical2BC* kmer_factory, uint64_t min_abundance, bool only_canonical, std::string output_path);

        void copy_content_from_another_hash_table(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* new_hasher, PointerHashTableCanonical* old_hash_table);

};