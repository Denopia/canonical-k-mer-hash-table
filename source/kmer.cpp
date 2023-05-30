#include "kmer.hpp"

// Constructor
OneCharacterAndPointerKMerAtomicFlag::OneCharacterAndPointerKMerAtomicFlag()
{
    data = 0ULL;
    my_flag.clear(std::memory_order_release);
}

// Destructor
OneCharacterAndPointerKMerAtomicFlag::~OneCharacterAndPointerKMerAtomicFlag(){}

// Acquire flag function
void OneCharacterAndPointerKMerAtomicFlag::acquire_lock()
{
    while(my_flag.test_and_set(std::memory_order_acquire));
}
// Release flag function
void OneCharacterAndPointerKMerAtomicFlag::release_lock()
{
    my_flag.clear(std::memory_order_release);
}

// Previous k-mer functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_predecessor_slot(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = data >> (64-38);
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return data >> (64-38);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_predecessor_slot(uint64_t previous_kmer_slot, bool atomic)
{
    
    // 2**38-1 = 274,877,906,943
    // Discard pointer mask
    //uint64_t discard_pointer_mask = 67108863;
    // 64 - 38 = 26
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        // Discard old pointer
        data = (data & uint64_t(67108863));
        // Set new pointer
        data = (data | (previous_kmer_slot << 26));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        // Discard old pointer
        data = (data & uint64_t(67108863));
        // Set new pointer
        data = (data | (previous_kmer_slot << 26));
    }
}

// Full data block functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_data(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = data;
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return data;
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_data(uint64_t new_data, bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data = new_data;
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = new_data;
    }
}

// Count functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_count(bool atomic)
{
    // max count 16,383
    //uint64_t count_mask = 16383;
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data >> 12) & uint64_t(16383));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data >> 12) & uint64_t(16383));
    }
    //return ((data >> 12) & uint64_t(16383));
}

void OneCharacterAndPointerKMerAtomicFlag::set_count(uint64_t new_count, bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        uint64_t discard_count_mask = ~(uint64_t(16383) << 12); 
        data = (data & discard_count_mask);
        data = (data | (new_count << 12));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        uint64_t discard_count_mask = ~(uint64_t(16383) << 12); 
        data = (data & discard_count_mask);
        data = (data | (new_count << 12));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::increase_count(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        uint64_t my_count = ((data >> 12) & uint64_t(16383));
        if (my_count != 16383)
        {
            // 2**12 = 4096 == 0001 0000 0000 0000
            data += 4096;
        }
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        uint64_t my_count = ((data >> 12) & uint64_t(16383));
        if (my_count != 16383)
        {
            // 2**12 = 4096 == 0001 0000 0000 0000
            data += 4096;
        }
    }
}

// Character functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_left_character(bool atomic)
{
    if(atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data>>10) & uint32_t(3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data>>10) & uint32_t(3));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_left_character(uint64_t new_character,bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //uint64_t discard_char_mask = ~((3) << 10); 
        // 1100 0000 0000 = 3072
        data = (data & (~(uint64_t(3072))));
        data = (data | (new_character << 10));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = (data & (~(uint64_t(3072))));
        data = (data | (new_character << 10));
    }
}

uint64_t OneCharacterAndPointerKMerAtomicFlag::get_right_character(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data>>8) & uint32_t(3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data>>8) & uint32_t(3));
    }
}

// DONE until this
void OneCharacterAndPointerKMerAtomicFlag::set_right_character(uint64_t new_character,bool atomic)
{ 
    // 0011 0000 0000 = 768
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data = (data & (~(uint64_t(768))));
        data = (data | (new_character << 8));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = (data & (~(uint64_t(768))));
        data = (data | (new_character << 8));
    }
}

// Getter data functions (1 bit flags)
bool OneCharacterAndPointerKMerAtomicFlag::is_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1)) == uint64_t(1));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1)) == uint64_t(1));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<1)) == uint64_t(1<<1));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<1)) == uint64_t(1<<1));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<2)) == uint64_t(1<<2));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<2)) == uint64_t(1<<2));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<3)) == uint64_t(1<<3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<3)) == uint64_t(1<<3));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<4)) == uint64_t(1<<4));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<4)) == uint64_t(1<<4));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<5)) == uint64_t(1<<5));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<5)) == uint64_t(1<<5));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<6)) == uint64_t(1<<6));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<6)) == uint64_t(1<<6));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<7)) == uint64_t(1<<7));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<7)) == uint64_t(1<<7));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_complete(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool leftnull = ((data & uint64_t(1<<2)) == uint64_t(1<<2));
        bool rightnull = ((data & uint64_t(1<<3)) == uint64_t(1<<3));
        my_flag.clear(std::memory_order_release);
        return ((!leftnull) && (!rightnull));
    }
    else
    {
        return ((!left_char_is_null(false)) && (!right_char_is_null(false)));
    }
}

// Setter data functions (1 bit flags)

void OneCharacterAndPointerKMerAtomicFlag::set_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<1);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<1);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 1);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 1);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<2);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<2);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 2);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 2);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<3);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<3);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 3);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 3);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<4);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<4);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 4);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 4);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<5);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<5);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 5);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 5);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<6);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<6);
    }
}
        
void OneCharacterAndPointerKMerAtomicFlag::unset_is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 6);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 6);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<7);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<7);
    }
}
        
void OneCharacterAndPointerKMerAtomicFlag::unset_is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 7);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 7);
    }
}


/*

// THIS IS FOR FLAGLESS K-MER WITH ATOMIC DATA VARIABLE

*/

//============
// Constructor
//============
OneCharacterAndPointerKMerAtomicVariable::OneCharacterAndPointerKMerAtomicVariable()
{
    data.store(0ULL, std::memory_order_release);
}

//===========
// Destructor
//===========
OneCharacterAndPointerKMerAtomicVariable::~OneCharacterAndPointerKMerAtomicVariable(){}


//==========================================
// Functions to fetch specific parts of data
//==========================================

uint64_t OneCharacterAndPointerKMerAtomicVariable::get_data()
{
    return data.load(std::memory_order_acquire);
}

uint64_t OneCharacterAndPointerKMerAtomicVariable::get_predecessor_slot()
{
    return data.load(std::memory_order_acquire) >> (64-38);
}

uint64_t OneCharacterAndPointerKMerAtomicVariable::get_count()
{
    return ((data.load(std::memory_order_acquire) >> 12) & uint64_t(16383));
}

uint64_t OneCharacterAndPointerKMerAtomicVariable::get_left_character()
{
    return ((data.load(std::memory_order_acquire) >> 10) & uint64_t(3));
}

uint64_t OneCharacterAndPointerKMerAtomicVariable::get_right_character()
{
    return ((data.load(std::memory_order_acquire) >> 8) & uint64_t(3));
}

bool OneCharacterAndPointerKMerAtomicVariable::is_occupied()
{
    return (data.load(std::memory_order_acquire) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::predecessor_exists()
{
    return ((data.load(std::memory_order_acquire) >> 1) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::left_char_is_null()
{
    return ((data.load(std::memory_order_acquire) >> 2) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::right_char_is_null()
{
    return ((data.load(std::memory_order_acquire) >> 3) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::canonical_during_insertion_self()
{
    return ((data.load(std::memory_order_acquire) >> 4) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::canonical_during_insertion_predecessor()
{
    return ((data.load(std::memory_order_acquire) >> 5) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::is_flagged_1()
{
    return ((data.load(std::memory_order_acquire) >> 6) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::is_flagged_2()
{
    return ((data.load(std::memory_order_acquire) >> 7) & uint64_t(1));
}

bool OneCharacterAndPointerKMerAtomicVariable::is_complete()
{
    return (left_char_is_null() && right_char_is_null());
}

//========================================
// Function to increase counter atomically
//========================================

void OneCharacterAndPointerKMerAtomicVariable::increase_count()
{
    // Get current data
    uint64_t current_data = data.load(std::memory_order_acquire);
    if (((current_data >> 12) & uint64_t(16383)) == uint64_t(16383))
        return;
    // Try to increase count
    //while(!data.compare_exchange_strong(current_data, kmod::modify_to_increase_count_by_one(current_data), std::memory_order_release, std::memory_order_relaxed))
    while(!data.compare_exchange_strong(current_data, kmod::modify_to_increase_count_by_one(current_data), std::memory_order_acq_rel, std::memory_order_relaxed))
    //while(!data.compare_exchange_weak(current_data, kmod::modify_to_increase_count_by_one(current_data), std::memory_order_release, std::memory_order_relaxed))
    {
        // If it was already at max, do nothing
        if (((current_data >> 12) & uint64_t(16383)) == uint64_t(16383))
            break;
    }
}




/*
// Constructor
OneCharacterAndPointerKMerAtomicVariable::OneCharacterAndPointerKMerAtomicVariable()
{
    data.store(0ULL);
}

// Destructor
OneCharacterAndPointerKMerAtomicVariable::~OneCharacterAndPointerKMerAtomicVariable(){}


// Previous k-mer functions
uint64_t OneCharacterAndPointerKMerAtomicVariable::get_predecessor_slot()
{
    return data.load() >> (64-38);
}

uint64_t OneCharacterAndPointerKMerAtomicVariable::set_predecessor_slot_to_data(uint64_t current_data, uint64_t previous_kmer_slot)
{
    // Discard old pointer
    current_data = (current_data & uint64_t(67108863));
    // Set new pointer
    current_data = (current_data | (previous_kmer_slot << 26));

    return current_data;
}

// Full data block functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_data()
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = data;
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return data;
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_data(uint64_t new_data, bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data = new_data;
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = new_data;
    }
}

// Count functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_count(bool atomic)
{
    // max count 16,383
    //uint64_t count_mask = 16383;
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data >> 12) & uint64_t(16383));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data >> 12) & uint64_t(16383));
    }
    //return ((data >> 12) & uint64_t(16383));
}

void OneCharacterAndPointerKMerAtomicFlag::set_count(uint64_t new_count, bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        uint64_t discard_count_mask = ~(uint64_t(16383) << 12); 
        data = (data & discard_count_mask);
        data = (data | (new_count << 12));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        uint64_t discard_count_mask = ~(uint64_t(16383) << 12); 
        data = (data & discard_count_mask);
        data = (data | (new_count << 12));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::increase_count(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        uint64_t my_count = ((data >> 12) & uint64_t(16383));
        if (my_count != 16383)
        {
            // 2**12 = 4096 == 0001 0000 0000 0000
            data += 4096;
        }
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        uint64_t my_count = ((data >> 12) & uint64_t(16383));
        if (my_count != 16383)
        {
            // 2**12 = 4096 == 0001 0000 0000 0000
            data += 4096;
        }
    }
}

// Character functions
uint64_t OneCharacterAndPointerKMerAtomicFlag::get_left_character(bool atomic)
{
    if(atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data>>10) & uint32_t(3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data>>10) & uint32_t(3));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_left_character(uint64_t new_character,bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //uint64_t discard_char_mask = ~((3) << 10); 
        // 1100 0000 0000 = 3072
        data = (data & (~(uint64_t(3072))));
        data = (data | (new_character << 10));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = (data & (~(uint64_t(3072))));
        data = (data | (new_character << 10));
    }
}

uint64_t OneCharacterAndPointerKMerAtomicFlag::get_right_character(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        uint64_t rv = ((data>>8) & uint32_t(3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data>>8) & uint32_t(3));
    }
}

// DONE until this
void OneCharacterAndPointerKMerAtomicFlag::set_right_character(uint64_t new_character,bool atomic)
{ 
    // 0011 0000 0000 = 768
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data = (data & (~(uint64_t(768))));
        data = (data | (new_character << 8));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data = (data & (~(uint64_t(768))));
        data = (data | (new_character << 8));
    }
}

// Getter data functions (1 bit flags)
bool OneCharacterAndPointerKMerAtomicFlag::is_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1)) == uint64_t(1));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1)) == uint64_t(1));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<1)) == uint64_t(1<<1));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<1)) == uint64_t(1<<1));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<2)) == uint64_t(1<<2));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<2)) == uint64_t(1<<2));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<3)) == uint64_t(1<<3));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<3)) == uint64_t(1<<3));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<4)) == uint64_t(1<<4));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<4)) == uint64_t(1<<4));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<5)) == uint64_t(1<<5));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<5)) == uint64_t(1<<5));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<6)) == uint64_t(1<<6));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<6)) == uint64_t(1<<6));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool rv = ((data & uint64_t(1<<7)) == uint64_t(1<<7));
        my_flag.clear(std::memory_order_release);
        return rv;
    }
    else
    {
        return ((data & uint64_t(1<<7)) == uint64_t(1<<7));
    }
}

bool OneCharacterAndPointerKMerAtomicFlag::is_complete(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        //while(my_flag.test(std::memory_order_acquire));
        bool leftnull = ((data & uint64_t(1<<2)) == uint64_t(1<<2));
        bool rightnull = ((data & uint64_t(1<<3)) == uint64_t(1<<3));
        my_flag.clear(std::memory_order_release);
        return ((!leftnull) && (!rightnull));
    }
    else
    {
        return ((!left_char_is_null(false)) && (!right_char_is_null(false)));
    }
}

// Setter data functions (1 bit flags)

void OneCharacterAndPointerKMerAtomicFlag::set_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_occupied(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1));
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1));
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<1);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<1);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_predecessor_exists(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 1);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 1);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<2);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<2);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_left_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 2);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 2);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<3);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<3);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_right_char_is_null(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 3);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 3);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<4);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<4);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_canonical_during_insertion_self(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 4);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 4);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<5);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<5);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::unset_canonical_during_insertion_predecessor(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 5);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 5);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<6);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<6);
    }
}
        
void OneCharacterAndPointerKMerAtomicFlag::unset_is_flagged_1(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 6);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 6);
    }
}

void OneCharacterAndPointerKMerAtomicFlag::set_is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data |= (uint64_t(1)<<7);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data |= (uint64_t(1)<<7);
    }
}
        
void OneCharacterAndPointerKMerAtomicFlag::unset_is_flagged_2(bool atomic)
{
    if (atomic)
    {
        while(my_flag.test_and_set(std::memory_order_acquire));
        data &= ~(uint64_t(1) << 7);
        my_flag.clear(std::memory_order_release);
    }
    else
    {
        data &= ~(uint64_t(1) << 7);
    }
}
*/