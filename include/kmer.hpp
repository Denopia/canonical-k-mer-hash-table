#include <cstdint>
#include <limits>
#include <iostream>

#pragma once


class OneCharacterAndPointerKMer
{
    private:
        uint64_t data;
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
    public:
        // Constructor
        OneCharacterAndPointerKMer();
        
        // Destructor
        ~OneCharacterAndPointerKMer();
        
        // Previous k-mer functions
        uint64_t get_predecessor_slot();
        void set_predecessor_slot(uint64_t predecessor_slot);
        
        // Data block functions
        uint64_t get_data();
        void set_data(uint64_t new_data);
        
        // Count functions
        uint64_t get_count();
        void set_count(uint64_t new_count);
        void increase_count();
        
        // Character functions
        uint64_t get_left_character();
        void set_left_character(uint64_t new_left_character);
        uint64_t get_right_character();
        void set_right_character(uint64_t new_right_character);
        
        // Getter data functions
        bool is_occupied();
        bool predecessor_exists();
        bool left_char_is_null();
        bool right_char_is_null();
        bool canonical_during_insertion_self();
        bool canonical_during_insertion_predecessor();
        bool is_flagged_1();
        bool is_flagged_2();
        bool is_complete(); // inferred from 4 and 8
        
        // Setter data functions (1 bit flags)
        void set_occupied();
        void unset_occupied();
        void set_predecessor_exists();
        void unset_predecessor_exists();
        void set_left_char_is_null();
        void unset_left_char_is_null();
        void set_right_char_is_null();
        void unset_right_char_is_null();
        void set_canonical_during_insertion_self();
        void unset_canonical_during_insertion_self();
        void set_canonical_during_insertion_predecessor();
        void unset_canonical_during_insertion_predecessor();
        void set_is_flagged_1();
        void unset_is_flagged_1();
        void set_is_flagged_2();
        void unset_is_flagged_2();
        
};