#include "kmer.hpp"


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