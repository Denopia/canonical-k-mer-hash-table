#include <cstdint>
#include <limits>
#include <iostream>
#include <atomic>
#include "functions_kmer_mod.hpp"

#pragma once

class BasicAtomicKMer
{
    public:
        std::atomic<uint64_t> kmer;
        std::atomic<uint64_t> count;
        BasicAtomicKMer(): kmer(0), count(0) {};
        ~BasicAtomicKMer() {};
};


class OneCharacterAndPointerKMerAtomicFlag
{
    
    public:
        //std::atomic<uint64_t> data;
        uint64_t data;
        std::atomic_flag my_flag;
        /*
            Data (from right to left):
            * 38 bits for pointer (MAX 274,877,906,944 pointers: if at max -> uses 2 terabytes)
            * 14 bits for count (MAX 16,383)
            * 2 bits for left character
            * 2 bits for right character
            * 8 bits for various flags (from right to left)
                - (128) 1 bit for free flag 2 (marks anything that is needed)       
                -  (64) 1 bit for free flag 1 (marks anything that is needed)
                -  (32) 1 bit for predecessor was canonical in read when inserted
                -  (16) 1 bit for self was canonical in read when inserted
                -   (8) 1 bit for right character is null
                -   (4) 1 bit for left character is null
                -   (2) 1 bit for predecessor k-mer exists (1 = exists, 0 = does not exist)
                -   (1) 1 bit for occupied (1 = occupied, 0 = free)
        */
        // Constructor
        OneCharacterAndPointerKMerAtomicFlag();
        
        // Destructor
        ~OneCharacterAndPointerKMerAtomicFlag();
        
        // Flag functions
        void acquire_lock();
        void release_lock();

        // Previous k-mer functions
        uint64_t get_predecessor_slot(bool atomic = true);
        void set_predecessor_slot(uint64_t predecessor_slot, bool atomic = true);
        
        // Data block functions
        uint64_t get_data(bool atomic = true);
        void set_data(uint64_t new_data,bool atomic = true);
        
        // Count functions
        uint64_t get_count(bool atomic = true);
        void set_count(uint64_t new_count,bool atomic = true);
        void increase_count(bool atomic = true);
        
        // Character functions
        uint64_t get_left_character(bool atomic = true);
        void set_left_character(uint64_t new_left_character,bool atomic = true);
        uint64_t get_right_character(bool atomic = true);
        void set_right_character(uint64_t new_right_character,bool atomic = true);
        
        // Getter data functions
        bool is_occupied(bool atomic = true);
        bool predecessor_exists(bool atomic = true);
        bool left_char_is_null(bool atomic = true);
        bool right_char_is_null(bool atomic = true);
        bool canonical_during_insertion_self(bool atomic = true);
        bool canonical_during_insertion_predecessor(bool atomic = true);
        bool is_flagged_1(bool atomic = true);
        bool is_flagged_2(bool atomic = true);
        bool is_complete(bool atomic = true); // inferred from 4 and 8
        
        // Setter data functions (1 bit flags)
        void set_occupied(bool atomic = true);
        void unset_occupied(bool atomic = true);
        void set_predecessor_exists(bool atomic = true);
        void unset_predecessor_exists(bool atomic = true);
        void set_left_char_is_null(bool atomic = true);
        void unset_left_char_is_null(bool atomic = true);
        void set_right_char_is_null(bool atomic = true);
        void unset_right_char_is_null(bool atomic = true);
        void set_canonical_during_insertion_self(bool atomic = true);
        void unset_canonical_during_insertion_self(bool atomic = true);
        void set_canonical_during_insertion_predecessor(bool atomic = true);
        void unset_canonical_during_insertion_predecessor(bool atomic = true);
        void set_is_flagged_1(bool atomic = true);
        void unset_is_flagged_1(bool atomic = true);
        void set_is_flagged_2(bool atomic = true);
        void unset_is_flagged_2(bool atomic = true);
        
};


class OneCharacterAndPointerKMerAtomicVariable
{
    
    public:
        std::atomic<uint64_t> data;
        /*
            Data (from right to left):
            * 38 bits for pointer (MAX 274,877,906,944 pointers: if at max -> uses 2 terabytes)
            * 14 bits for count (MAX 16,383)
            * 2 bits for left character
            * 2 bits for right character
            * 8 bits for various flags (from right to left)
                - (128) 1 bit for free flag 2 (marks anything that is needed)       
                -  (64) 1 bit for free flag 1 (marks anything that is needed)
                -  (32) 1 bit for predecessor was canonical in read when inserted
                -  (16) 1 bit for self was canonical in read when inserted
                -   (8) 1 bit for right character is null
                -   (4) 1 bit for left character is null
                -   (2) 1 bit for predecessor k-mer exists (1 = exists, 0 = does not exist)
                -   (1) 1 bit for occupied (1 = occupied, 0 = free)
        */
        // Constructor
        OneCharacterAndPointerKMerAtomicVariable();
        
        // Destructor
        ~OneCharacterAndPointerKMerAtomicVariable();
        
        // Getters
        uint64_t get_data();
        uint64_t get_predecessor_slot();
        uint64_t get_count();
        uint64_t get_left_character();
        uint64_t get_right_character();
        bool is_occupied();
        bool predecessor_exists();
        bool left_char_is_null();
        bool right_char_is_null();
        bool canonical_during_insertion_self();
        bool canonical_during_insertion_predecessor();
        bool is_flagged_1();
        bool is_flagged_2();
        bool is_complete(); // inferred from 4 and 8
        
        // Counter increaser (+1)
        void increase_count();

};
