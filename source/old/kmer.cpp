#include "kmer.hpp"

//##################### THE NEW COMPACT VERSION STARTS #######################################
// Constructor
OneCharacterAndPointerKMer::OneCharacterAndPointerKMer(): data(0){}

// Destructor
OneCharacterAndPointerKMer::~OneCharacterAndPointerKMer(){}

// Previous k-mer functions
uint64_t OneCharacterAndPointerKMer::get_previous_kmer_slot()
{
    return data >> (64-38);
    //return previous_kmer_slot;
}

void OneCharacterAndPointerKMer::set_previous_kmer_slot(uint64_t previous_kmer_slot)
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

uint64_t OneCharacterAndPointerKMer::get_character(bool take_left_char, bool flip_char)
{
    if (take_left_char){
        if (flip_char){
            return uint64_t(3)-get_left_character();
        } else {
            return get_left_character();
        }
    } else {
        if (flip_char){
            return uint64_t(3)-get_right_character();
        } else {
            return get_right_character();
        }
    }
}

uint64_t OneCharacterAndPointerKMer::get_canonical_right_character()
{
    return get_right_character();
}

// Getter data functions (1 bit flags)
bool OneCharacterAndPointerKMer::is_occupied()
{
    return ((data & uint64_t(1)) == uint64_t(1));
    //return ((data & uint64_t(1)) == uint64_t(1));
}

bool OneCharacterAndPointerKMer::is_complete()
{
    return ((data & uint64_t(1<<1)) == uint64_t(1<<1));
    //return ((data & uint64_t(2)) == uint64_t(2)); 
}

bool OneCharacterAndPointerKMer::prev_kmer_exists()
{
    return ((data & uint64_t(1<<2)) == uint64_t(1<<2));
    //return ((data & uint64_t(4)) == uint64_t(4));
}

bool OneCharacterAndPointerKMer::check_starts_from_left()
{
    return ((data & uint64_t(1<<3)) == uint64_t(1<<3));
    //return ((data & uint64_t(8)) == uint64_t(8));
}

bool OneCharacterAndPointerKMer::same_orientation_as_previous()
{
    return ((data & uint64_t(1<<4)) == uint64_t(1<<4));
    //return ((data & uint64_t(16)) == uint64_t(16));
}

bool OneCharacterAndPointerKMer::is_in_main_array()
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



// Setter data functions (1 bit flags)

void OneCharacterAndPointerKMer::set_occupied()
{
    data |= (uint64_t(1));
}

void OneCharacterAndPointerKMer::unset_occupied()
{
    data &= ~(uint64_t(1));
}

void OneCharacterAndPointerKMer::set_complete()
{
    data |= (uint64_t(1)<<1);
}

void OneCharacterAndPointerKMer::unset_complete()
{
    data &= ~(uint64_t(1) << 1);
}

void OneCharacterAndPointerKMer::set_prev_kmer_exists()
{
    data |= (uint64_t(1)<<2);
}

void OneCharacterAndPointerKMer::unset_prev_kmer_exists()
{
    data &= ~(uint64_t(1) << 2);
}

void OneCharacterAndPointerKMer::set_check_starts_from_left()
{
    data |= (uint64_t(1)<<3);
}

void OneCharacterAndPointerKMer::unset_check_starts_from_left()
{
    data &= ~(uint64_t(1) << 3);
}

void OneCharacterAndPointerKMer::set_same_orientation_as_previous()
{
    data |= (uint64_t(1)<<4);
}

void OneCharacterAndPointerKMer::unset_same_orientation_as_previous()
{
    data &= ~(uint64_t(1) << 4);
}

void OneCharacterAndPointerKMer::set_in_main_array()
{
    data |= (uint64_t(1)<<5);
}

void OneCharacterAndPointerKMer::unset_in_main_array()
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
//##################### THE NEW COMPACT VERSION ENDS #######################################


//##################### THE ORIGINAL WORKING VERSION STARTS #######################################
/*
// Constructor
OneCharacterAndPointerKMer::OneCharacterAndPointerKMer(): previous_kmer_slot(0), data(0), count(0){}

// Destructor
OneCharacterAndPointerKMer::~OneCharacterAndPointerKMer(){}

// Previous k-mer functions

uint32_t OneCharacterAndPointerKMer::get_previous_kmer_slot()
{
    return previous_kmer_slot;
}

void OneCharacterAndPointerKMer::set_previous_kmer_slot(uint32_t previous)
{
    previous_kmer_slot = previous;
}

// Count functions

uint32_t OneCharacterAndPointerKMer::get_count()
{
    return count;
}

uint32_t OneCharacterAndPointerKMer::get_data()
{
    return data;
}

void OneCharacterAndPointerKMer::set_data(uint32_t new_data)
{
    data = new_data;
}

void OneCharacterAndPointerKMer::set_count(uint32_t new_count)
{
    count = new_count;
}

void OneCharacterAndPointerKMer::increase_count()
{
    if (is_at_max_count())
    {
        std::cout << "k-mer count already at MAX\n";
        return;
    }

    count += 1;

    if (count == std::numeric_limits<uint32_t>::max())
        set_at_max_count();
}

// Character functions

uint32_t OneCharacterAndPointerKMer::get_character()
{
    return data & uint32_t(3);
}

void OneCharacterAndPointerKMer::set_character(uint32_t new_character)
{
    data &= (~uint32_t(3));
    data |= (new_character&uint32_t(3));
}

// Getter data functions (1 bit flags)

bool OneCharacterAndPointerKMer::is_occupied()
{
    return ((data & uint32_t(1<<2)) == uint32_t(1<<2));
    //return ((data & uint32_t(4)) == uint32_t(4));
}

bool OneCharacterAndPointerKMer::is_complete()
{
    return ((data & uint32_t(1<<3)) == uint32_t(1<<3));
    //return ((data & uint32_t(8)) == uint32_t(8)); 
}

bool OneCharacterAndPointerKMer::prev_kmer_exists()
{
    return ((data & uint32_t(1<<4)) == uint32_t(1<<4));
    //return ((data & uint32_t(16)) == uint32_t(16));
}

bool OneCharacterAndPointerKMer::is_from_forward_strand()
{
    return ((data & uint32_t(1<<5)) == uint32_t(1<<5));
    //return ((data & uint32_t(32)) == uint32_t(32));
}

bool OneCharacterAndPointerKMer::is_at_max_count()
{
    return ((data & uint32_t(1<<6)) == uint32_t(1<<6));
    //return ((data & uint32_t(64)) == uint32_t(64));
}

bool OneCharacterAndPointerKMer::is_in_main_array()
{
    return ((data & uint32_t(1<<7)) == uint32_t(1<<7));
    //return ((data & uint32_t(128)) == uint32_t(128));
}

bool OneCharacterAndPointerKMer::is_written_in_output()
{
    return ((data & uint32_t(1<<8)) == uint32_t(1<<8));
    //return ((data & uint32_t(256)) == uint32_t(256));
}

bool OneCharacterAndPointerKMer::is_flagged_1()
{
    return ((data & uint32_t(1<<9)) == uint32_t(1<<9));
    //return ((data & uint32_t(512)) == uint32_t(256));
}

bool OneCharacterAndPointerKMer::is_flagged_2()
{
    return ((data & uint32_t(1<<10)) == uint32_t(1<<10));
}

// Setter data functions (1 bit flags)

void OneCharacterAndPointerKMer::set_occupied()
{
    data |= (uint32_t(1)<<2);
}

void OneCharacterAndPointerKMer::unset_occupied()
{
    data &= ~(uint32_t(1) << 2);
}

void OneCharacterAndPointerKMer::set_complete()
{
    data |= (uint32_t(1)<<3);
}

void OneCharacterAndPointerKMer::unset_complete()
{
    data &= ~(uint32_t(1) << 3);
}

void OneCharacterAndPointerKMer::set_prev_kmer_exists()
{
    data |= (uint32_t(1)<<4);
}

void OneCharacterAndPointerKMer::unset_prev_kmer_exists()
{
    data &= ~(uint32_t(1) << 4);
}

void OneCharacterAndPointerKMer::set_from_forward_strand()
{
    data |= (uint32_t(1)<<5);
}

void OneCharacterAndPointerKMer::unset_from_forward_strand()
{
    data &= ~(uint32_t(1) << 5);
}

void OneCharacterAndPointerKMer::set_at_max_count()
{
    data |= (uint32_t(1)<<6);
}

void OneCharacterAndPointerKMer::unset_at_max_count()
{
    data &= ~(uint32_t(1) << 6);
}

void OneCharacterAndPointerKMer::set_in_main_array()
{
    data |= (uint32_t(1)<<7);
}

void OneCharacterAndPointerKMer::unset_in_main_array()
{
    data &= ~(uint32_t(1) << 7);
}

void OneCharacterAndPointerKMer::set_is_written_in_output()
{
    data |= (uint32_t(1)<<8);
}

void OneCharacterAndPointerKMer::unset_is_written_in_output()
{
    data &= ~(uint32_t(1) << 8);
}

void OneCharacterAndPointerKMer::set_is_flagged_1()
{
    data |= (uint32_t(1)<<9);
}
        
void OneCharacterAndPointerKMer::unset_is_flagged_1()
{
    data &= ~(uint32_t(1) << 9);
}

void OneCharacterAndPointerKMer::set_is_flagged_2()
{
    data |= (uint32_t(1)<<10);
}
        
void OneCharacterAndPointerKMer::unset_is_flagged_2()
{
    data &= ~(uint32_t(1) << 10);
}
*/
//##################### THE ORIGINAL WORKING VERSION ENDS #######################################