#include <vector>
#include <math.h>
#include <cstdint>
#include <iostream>
//#include "kmer.hpp"
#include "functions_strings.hpp"

#pragma once

// 2 bits per character version of k-mer factory
class KMerFactory2BC 
{
    public:

        uint64_t* blocks; // V
        uint64_t used_left_block_mask; // V

        uint64_t right_block_right_char_mask; // V
        uint64_t left_block_left_char_mask; // V

        uint64_t pushed_off_character; // V

        uint64_t right_char_mask; // V
        uint64_t left_char_mask; // V

        int character_bits;
        int kmer_length; // V
        int characters_stored; // V
        
        int bits_in_last_block; // V
        int number_of_blocks; // V


        // Constructor
        KMerFactory2BC(uint64_t k);

        // Destructor
        ~KMerFactory2BC();

        // Resets the factory
        void reset();

        // Return the number of stored characters
        int get_number_of_stored_characters();

        // Check if the number of stored characters is equal to k
        bool current_kmer_is_real();

        // Push new character to the factory
        void push_new_character(char c);

        // Push new character to the factory in integer format
        void push_new_integer(uint64_t c);

        // Push a new full block of characters
        void push_new_block(uint64_t b);

        // Return the leftmost character of the leftmost block
        // This is also the oldest character
        uint64_t get_leftmost_character();

        // Return the rightmost character of the rightmost block
        // This is also the newest character
        uint64_t get_rightmost_character();

        uint64_t get_newest_character();

        // Returns the character that was most recently pushed off from the leftmost block
        // aka the character on the left from the leftmost character
        uint64_t get_pushed_off_character();

};

// 2 bits per character version of k-mer factory with canonical functionality
class KMerFactoryCanonical2BC 
{
    public:

        int character_bits;
        int kmer_length; // V
        int characters_stored; // V
        uint64_t pushed_off_character_forward; // V

        uint64_t right_char_mask; // V
        uint64_t left_char_mask; // V

        int bits_in_last_block; // V
        int number_of_blocks; // V

        uint64_t right_block_right_char_mask; // V
        uint64_t left_block_left_char_mask; // V

        uint64_t used_left_block_mask; // V
        
        uint64_t* blocks_forward; // V
        uint64_t* blocks_backward; // V

        bool forward_is_canonical;
        bool previous_forward_was_canonical;
        bool previous_kmer_exists;
        

        // Constructor
        KMerFactoryCanonical2BC(uint64_t k);

        // Destructor
        ~KMerFactoryCanonical2BC();

        // Resets the factory
        void reset();

        // Find out if the forward version of the current k-mer is canonical
        bool forward_kmer_is_canonical();

        // Return the number of stored characters
        int get_number_of_stored_characters();

        // Check if the number of stored characters is equal to k
        bool current_kmer_is_real();

        // Push new character to the factory
        void push_new_character(char c);

        // Push new character to the factory
        void push_new_integer(uint64_t c);

        // Return the leftmost character of the leftmost block
        // This is also the oldest character
        uint64_t get_forward_leftmost_character();

        // Return the rightmost character of the rightmost block
        // This is also the newest character
        uint64_t get_forward_rightmost_character();

        // Returns the character that was most recently pushed off from the leftmost block
        // aka the character on the left from the leftmost character
        uint64_t get_forward_pushed_off_character();

        uint64_t get_forward_newest_character();

        uint64_t get_canonical_block(uint64_t i);

        uint64_t get_noncanonical_block(uint64_t i);

        uint64_t get_canonical_right_character();

        uint64_t get_canonical_left_character();

};
