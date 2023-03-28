#include "kmer.hpp"

// Constructor
OneCharacterAndPointerKMer::OneCharacterAndPointerKMer(): data(0){}

// Destructor
OneCharacterAndPointerKMer::~OneCharacterAndPointerKMer(){}

// Previous k-mer functions
uint64_t OneCharacterAndPointerKMer::get_predecessor_slot()
{
    return data >> (64-38);
}

void OneCharacterAndPointerKMer::set_predecessor_slot(uint64_t previous_kmer_slot)
{
    // 2**38-1 = 274,877,906,943
    // Discard pointer mask
    //uint64_t discard_pointer_mask = 67108863;
    // 64 - 38 = 26 
    // Discard old pointer
    data = (data & uint64_t(67108863));
    // Set new pointer
    data = (data | (previous_kmer_slot << 26));
}

// Full data block functions
uint64_t OneCharacterAndPointerKMer::get_data()
{
    return data;
}

void OneCharacterAndPointerKMer::set_data(uint64_t new_data)
{
    data = new_data;
}

// Count functions
uint64_t OneCharacterAndPointerKMer::get_count()
{
    // max count 16,383
    //uint64_t count_mask = 16383;
    return ((data >> 12) & uint64_t(16383));
}

void OneCharacterAndPointerKMer::set_count(uint64_t new_count)
{
    uint64_t discard_count_mask = ~(uint64_t(16383) << 12); 
    data = (data & discard_count_mask);
    data = (data | (new_count << 12));
}

void OneCharacterAndPointerKMer::increase_count()
{
    if (get_count() != 16383)
    {
        // 2**12 = 4096 == 0001 0000 0000 0000
        data += 4096;
    }
}

// Character functions
uint64_t OneCharacterAndPointerKMer::get_left_character()
{
    return ((data>>10) & uint32_t(3));
}

void OneCharacterAndPointerKMer::set_left_character(uint64_t new_character)
{
    //uint64_t discard_char_mask = ~((3) << 10); 
    // 1100 0000 0000 = 3072
    data = (data & (~(uint64_t(3072))));
    data = (data | (new_character << 10));
}

uint64_t OneCharacterAndPointerKMer::get_right_character()
{
    return ((data>>8) & uint32_t(3));
}

void OneCharacterAndPointerKMer::set_right_character(uint64_t new_character)
{ 
    // 0011 0000 0000 = 768
    data = (data & (~(uint64_t(768))));
    data = (data | (new_character << 8));
}

// Getter data functions (1 bit flags)
bool OneCharacterAndPointerKMer::is_occupied()
{
    return ((data & uint64_t(1)) == uint64_t(1));
    //return ((data & uint64_t(1)) == uint64_t(1));
}

bool OneCharacterAndPointerKMer::predecessor_exists()
{
    return ((data & uint64_t(1<<1)) == uint64_t(1<<1));
}

bool OneCharacterAndPointerKMer::left_char_is_null()
{
    return ((data & uint64_t(1<<2)) == uint64_t(1<<2));
    //return ((data & uint64_t(4)) == uint64_t(4));
}

bool OneCharacterAndPointerKMer::right_char_is_null()
{
    return ((data & uint64_t(1<<3)) == uint64_t(1<<3));
    //return ((data & uint64_t(8)) == uint64_t(8));
}

bool OneCharacterAndPointerKMer::canonical_during_insertion_self()
{
    return ((data & uint64_t(1<<4)) == uint64_t(1<<4));
    //return ((data & uint64_t(16)) == uint64_t(16));
}

bool OneCharacterAndPointerKMer::canonical_during_insertion_predecessor()
{
    return ((data & uint64_t(1<<5)) == uint64_t(1<<5));
    //return ((data & uint64_t(32)) == uint64_t(32));
}

bool OneCharacterAndPointerKMer::is_flagged_1()
{
    return ((data & uint64_t(1<<6)) == uint64_t(1<<6));
    //return ((data & uint64_t(64)) == uint64_t(256));
}

bool OneCharacterAndPointerKMer::is_flagged_2()
{
    return ((data & uint64_t(1<<7)) == uint64_t(1<<7));
    //return ((data & uint64_t(128)) == uint64_t(128));
}

bool OneCharacterAndPointerKMer::is_complete()
{
    return ((!left_char_is_null()) && (!right_char_is_null()));
}

// Setter data functions (1 bit flags)

void OneCharacterAndPointerKMer::set_occupied()
{
    data |= (uint64_t(1));
}

void OneCharacterAndPointerKMer::unset_occupied()
{
    data &= ~(uint64_t(1));
}

void OneCharacterAndPointerKMer::set_predecessor_exists()
{
    data |= (uint64_t(1)<<1);
}

void OneCharacterAndPointerKMer::unset_predecessor_exists()
{
    data &= ~(uint64_t(1) << 1);
}

void OneCharacterAndPointerKMer::set_left_char_is_null()
{
    data |= (uint64_t(1)<<2);
}

void OneCharacterAndPointerKMer::unset_left_char_is_null()
{
    data &= ~(uint64_t(1) << 2);
}

void OneCharacterAndPointerKMer::set_right_char_is_null()
{
    data |= (uint64_t(1)<<3);
}

void OneCharacterAndPointerKMer::unset_right_char_is_null()
{
    data &= ~(uint64_t(1) << 3);
}

void OneCharacterAndPointerKMer::set_canonical_during_insertion_self()
{
    data |= (uint64_t(1)<<4);
}

void OneCharacterAndPointerKMer::unset_canonical_during_insertion_self()
{
    data &= ~(uint64_t(1) << 4);
}

void OneCharacterAndPointerKMer::set_canonical_during_insertion_predecessor()
{
    data |= (uint64_t(1)<<5);
}

void OneCharacterAndPointerKMer::unset_canonical_during_insertion_predecessor()
{
    data &= ~(uint64_t(1) << 5);
}

void OneCharacterAndPointerKMer::set_is_flagged_1()
{
    data |= (uint64_t(1)<<6);
}
        
void OneCharacterAndPointerKMer::unset_is_flagged_1()
{
    data &= ~(uint64_t(1) << 6);
}

void OneCharacterAndPointerKMer::set_is_flagged_2()
{
    data |= (uint64_t(1)<<7);
}
        
void OneCharacterAndPointerKMer::unset_is_flagged_2()
{
    data &= ~(uint64_t(1) << 7);
}
