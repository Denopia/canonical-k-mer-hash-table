#include <vector>
#include <cmath>
#include <cstdint>
#include <iostream>
#include "functions_strings.hpp"

#pragma once

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
        explicit KMerFactoryCanonical2BC(int k);

        // Destructor
        ~KMerFactoryCanonical2BC();

    //TODO you can inline the functions below by moving their
    // definitions here and inserting the "inline" annotation.
    // Inlining should improve the performance

        // Check if the number of stored characters is equal to k
        [[nodiscard]] bool current_kmer_is_real() const;

        // Find out if the forward version of the current k-mer is canonical
        [[nodiscard]] bool forward_kmer_is_canonical() const;

        // Find out if the forward version of the previous k-mer was canonical
        [[nodiscard]] bool previous_forward_kmer_was_canonical() const;

        [[nodiscard]] bool previous_kmer_existed() const;

        // Return the number of stored characters
        [[nodiscard]] int get_number_of_stored_characters() const;

        // Resets the factory
        void reset();

        // Push new character to the factory
        void push_new_character(char c);

        // Push new character to the factory
        void push_new_integer(uint64_t c);

        // Return the leftmost character of the leftmost block
        // This is also the oldest character
        [[nodiscard]] uint64_t get_forward_leftmost_character() const;

        // Return the rightmost character of the rightmost block
        // This is also the newest character
        [[nodiscard]] uint64_t get_forward_rightmost_character() const;

        // Returns the character that was most recently pushed off from the leftmost block
        // aka the character on the left from the leftmost character
        [[nodiscard]] uint64_t get_forward_pushed_off_character() const;

        [[nodiscard]] uint64_t get_forward_newest_character() const;

        [[nodiscard]] uint64_t get_forward_block(uint64_t i) const;

        [[nodiscard]] uint64_t get_backward_block(uint64_t i) const;

        [[nodiscard]] uint64_t get_rightmost_forward_block() const;

        [[nodiscard]] uint64_t get_rightmost_backward_block() const;

        [[nodiscard]] uint64_t get_canonical_block(uint64_t i) const;

        [[nodiscard]] uint64_t get_noncanonical_block(uint64_t i) const;

        [[nodiscard]] uint64_t get_forward_char_at_position(int i) const;

        [[nodiscard]] uint64_t get_backward_char_at_position(int i) const;
};
