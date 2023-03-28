#include <cstdint>
#include <limits>
#include <iostream>


class OneCharacterAndPointerKMer
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
        OneCharacterAndPointerKMer();
        
        // Destructor
        ~OneCharacterAndPointerKMer();

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