#include <cstdint>
#include <limits>
#include <iostream>


class OneCharacterAndPointerKMer
{
    private:

        uint64_t data;

        // *** OLD ***
        /*
            Data (from right to left):
            * 39 bits for pointer (MAX 549,755,813,888, if at max -> uses 4 terabytes)
            * 14 bits for count (MAX 16,384)
            * 2 bits for character
            * 9 bits for various flags (from right to left)
                - (256) 1 bit for free flag 2 (marks anything that is needed)
                - (128) 1 bit for free flag 1 (marks anything that is needed)
                -  (64) 1 bit for written in output (1 = already written, 0 = not written yet)
                -  (32) 1 bit for in main array/in temp array (1 = in main array, 0 = in temporary array)
                -  (16) 1 bit for max count (1 = at max count, 0 = not at max count)
                -   (8) 1 bit for forward strand (1 = forward, 0 = reverse)
                -   (4) 1 bit for previous k-mer exists (1 = exists, 0 = does not exist)
                -   (2) 1 bit for complete/partial k-mer (1 = complete, 0 = partial)
                -   (1) 1 bit for occupied (1 = occupied, 0 = free)
        */


        // *** NEW ***
        /*
            Data (from right to left):
            * 38 bits for pointer (MAX 274,877,906,944 pointers: if at max -> uses 2 terabytes)
            * 14 bits for count (MAX 16,383)
            * 2 bits for left character
            * 2 bits for right character
            * 8 bits for various flags (from right to left)
                - (128) 1 bit for free flag 2 (marks anything that is needed)       
                -  (64) 1 bit for free flag 1 (marks anything that is needed)
                -  (32) 1 bit for in main array/in temp array (1 = in main array, 0 = in temporary array)
                -  (16) 1 bit for marking orientation w.r.t. the previous k-mer (1 = same, 0 = different)
                -   (8) 1 bit for marking from which end the check for this k-mer should be started (1 = from left, 0 = from right)
                -   (4) 1 bit for previous k-mer exists (1 = exists, 0 = does not exist)
                -   (2) 1 bit for complete/partial k-mer (1 = complete, 0 = partial)
                -   (1) 1 bit for occupied (1 = occupied, 0 = free)
        */

    public:

        // Constructor
        OneCharacterAndPointerKMer();
        
        // Destructor
        ~OneCharacterAndPointerKMer();

        // Previous k-mer functions

        uint64_t get_previous_kmer_slot();

        void set_previous_kmer_slot(uint64_t previous_kmer_slot);

        // Count functions

        uint64_t get_count();

        uint64_t get_data();

        void set_count(uint64_t new_count);

        void increase_count();

        void set_data(uint64_t new_data);

        // Character functions

        uint64_t get_left_character();

        void set_left_character(uint64_t new_character);

        uint64_t get_right_character();

        void set_right_character(uint64_t new_character);

        uint64_t get_character(bool take_left_char, bool flip_char);

        uint64_t get_canonical_right_character();

        // Getter data functions (1 bit flags)

        bool is_occupied();

        bool is_complete();

        bool prev_kmer_exists();

        bool check_starts_from_left();

        bool same_orientation_as_previous();

        bool is_in_main_array();

        bool is_flagged_1();

        bool is_flagged_2();
        

        // Setter data functions (1 bit flags)

        void set_occupied();

        void unset_occupied();

        void set_complete();

        void unset_complete();

        void set_prev_kmer_exists();

        void unset_prev_kmer_exists();

        void set_check_starts_from_left();

        void unset_check_starts_from_left();

        void set_same_orientation_as_previous();

        void unset_same_orientation_as_previous();

        void set_in_main_array();

        void unset_in_main_array();

        void set_is_flagged_1();
        
        void unset_is_flagged_1();

        void set_is_flagged_2();
        
        void unset_is_flagged_2();
};


class OneCharacterAndPointerKMer_OLD
{
    private:
        uint32_t previous_kmer_slot;
        uint32_t count;
        uint32_t data;
        /*
            Data (from right to left):
            - 2 bits for character
            - 1 bit for occupied (1 = occupied, 0 = free)
            - 1 bit for complete/partial k-mer (1 = complete, 0 = partial)
            - 1 bit for previous k-mer exists (1 = exists, 0 = does not exist)
            - 1 bit for forward strand (1 = forward, 0 = reverse)
            - 1 bit for max count (1 = at max count, 0 = not at max count)
            - 1 bit for in main array/in temp array (1 = in main array, 0 = in temporary array)
            - 1 bit for written in output (1 = already written, 0 = not written yet)
            - 1 bit for free flag 1 (marks anything that is needed)
            - 1 bit for free flag 2 (marks anything that is needed)
        */

    public:

        // Constructor
        OneCharacterAndPointerKMer_OLD();
        
        // Destructor
        ~OneCharacterAndPointerKMer_OLD();

        // Previous k-mer functions

        uint32_t get_previous_kmer_slot();

        void set_previous_kmer_slot(uint32_t previous);

        // Count functions

        uint32_t get_count();

        uint32_t get_data();

        void set_count(uint32_t new_count);

        void increase_count();

        // Character functions

        uint32_t get_character();

        void set_character(uint32_t new_character);

        void set_data(uint32_t new_data);

        // Getter data functions (1 bit flags)

        bool is_occupied();

        bool is_complete();

        bool prev_kmer_exists();

        bool is_from_forward_strand();

        bool is_at_max_count();

        bool is_in_main_array();

        bool is_written_in_output();

        bool is_flagged_1();

        bool is_flagged_2();
        

        // Setter data functions (1 bit flags)

        void set_occupied();

        void unset_occupied();

        void set_complete();

        void unset_complete();

        void set_prev_kmer_exists();

        void unset_prev_kmer_exists();

        void set_from_forward_strand();

        void unset_from_forward_strand();

        void set_at_max_count();

        void unset_at_max_count();

        void set_in_main_array();

        void unset_in_main_array();
        
        void set_is_written_in_output();

        void unset_is_written_in_output();

        void set_is_flagged_1();
        
        void unset_is_flagged_1();

        void set_is_flagged_2();
        
        void unset_is_flagged_2();
};