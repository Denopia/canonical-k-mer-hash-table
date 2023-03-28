#include "kmer_hash_table.hpp"

// ==============================================================================================================


// === For pointer hash table 2 =========================================================================================


PointerHashTable2::PointerHashTable2(uint64_t s, uint64_t k, uint64_t b)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMer[size];
    bits_per_char = 2;
    inserted_items = 0;
    inserted_complete_kmers = 0;
    solid_kmers = 0;
    kmer_blocks = b;
    probe_hasher = new ProbeHasher1();
    probing_prime = mathfunctions::next_prime(uint64_t(std::floor(size/13.0)));
    //max_temp_slots = std::max(std::ceil(size/200),100.0);
    max_temp_slots = 100;
    touched_temp_slots = 0;
    temp_slots_in_use = 0;
    max_temp_slot_in_use = 0;
    smallest_unused_temp_slot = 0;
    temp_array = std::vector<uint64_t>(b*max_temp_slots, uint64_t(0));
    temp_free_slots = std::vector<uint8_t>(max_temp_slots, 1);
    //temp_free_slots = sdsl::bit_vector(max_temp_slots, 1);
}

PointerHashTable2::~PointerHashTable2()
{
    delete[] hash_table_array;
    delete probe_hasher;
}

uint64_t PointerHashTable2::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(); 
}

bool PointerHashTable2::kmer_in_slot_is_complete(uint64_t slot)
{
    return hash_table_array[slot].is_complete(); 
}

void PointerHashTable2::resize()
{
    std::cout << "Hash table resizing not implemented yet...\n";
    exit(1);
}

uint64_t PointerHashTable2::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t PointerHashTable2::get_number_of_inserted_complete_kmers()
{
    return inserted_complete_kmers;
}

uint64_t PointerHashTable2::get_number_of_solid_kmers()
{
    return solid_kmers;
}

// Check if the currently loaded k-mer is in the given slot
// Return -1 if not in slot
// Return 0 if in slot but it is in the temp array from the very beginning
// Return 1 if in slot
int PointerHashTable2::kmer_slot_check(KMerFactory2BC* kmer_factory, uint64_t slot)
{
    //std::cout << "CHECKING MAIN ARRAY SLOT " << slot << "\n";
    // If slot is empty, k-mer is not here
    if (!hash_table_array[slot].is_occupied())
    {
        //std::cout << "Was not occupied\n";
        return -1;
    }
    
    // If slot is occupied by a k-mer in temp array, return false // DOES NOT WORK CORRECTLY IF TWO SEQUENCES START WITH THE EXACT SAME K-MER, FIX LATER
    // Should be fixed now
    if (!hash_table_array[slot].is_in_main_array())
    {
        //std::cout << "Found direct link to temp \n";
        //return false;
        for (int i = 0; i < kmer_factory->number_of_blocks; i++)
        {
            if (kmer_factory->blocks[i] != temp_array[hash_table_array[slot].get_previous_kmer_slot()*kmer_factory->number_of_blocks + i])
            {
                return -1;
            }
        }
        //std::cout << "Found using FULL TEMP\n";
        return 0;
    }
        
    // If we pass the initial checks, start comparing characters one by one
    uint64_t dbi = kmer_factory->number_of_blocks-1; // which data block we are looking at
    uint64_t db = kmer_factory->blocks[dbi]; // store current data block content here
    uint64_t cbc = 0; // number of checked data block characters
    uint64_t mc = 0; // matched characters
    bool check_rest_from_temp = false;
    //uint64_t temp_check_kmer_factory_start = 0;
    uint64_t temp_array_slot = 0;
    while (mc < kmer_factory->get_number_of_stored_characters())
    {
        // If the current hash table slot is unoccupied the k-mer does not exist
        if (!hash_table_array[slot].is_occupied())
        {
            return -1;
        }
        if (!hash_table_array[slot].is_in_main_array())
        {
            //std::cout << "Found link to temp mid search, start hybrid search\n";
            check_rest_from_temp = true;
            temp_array_slot = hash_table_array[slot].get_previous_kmer_slot();
            //std::cout << "Hybrid check starts in slot " << temp_array_slot << "\n";
            break;
        }
        // Check that the characters match
        if (uint64_t(hash_table_array[slot].get_character()) != (db&uint64_t(3)))
        {
            //std::cout << "k-mer not in this slot, normal check\n";
            return -1;
        // Update to look at next characters if the current ones match
        } 
        else 
        {
            mc+=1;
            if (mc == kmer_factory->get_number_of_stored_characters())
                return 1;
            if (hash_table_array[slot].prev_kmer_exists()){
                slot = hash_table_array[slot].get_previous_kmer_slot();
            } else {
                //std::cout << "VERY SUS ERROR!!!!!!!!!!!!!!!!!!!!!\n";
                return -1;
            }
            cbc += 1;
            db >>= 2;
            if (cbc >= 32){
                cbc = 0;
                dbi -= 1;
                db = kmer_factory->blocks[dbi];
            }
        }
    }
    // If the rest of he characters must be checked from temp
    if (check_rest_from_temp)
    {

        uint64_t current_temp_array_position = temp_array_slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1;
        uint64_t temp_array_block = temp_array[current_temp_array_position];
        int unchecked_temp_block_bits = 64;
        //std::cout << "Hybrid check adjusted slot is " << current_temp_array_position << "\n";
        //std::cout << "and the block in question is " << temp_array_block << "\n";
        //std::cout << "and the kmer factory block is now " << db << "\n";
        while (mc < kmer_factory->get_number_of_stored_characters())
        {
            //std::cout << "COMPARING FACTORY " << (db&uint64_t(3)) << " WITH TEMP " << (temp_array_block&uint64_t(3)) << "\n";
            if ((db&uint64_t(3)) != (temp_array_block&uint64_t(3)))
            {
                //std::cout << "k-mer not in this slot, hybrid check\n";
                return -1;
            }
            else
            {
                mc+=1;
                if (mc == kmer_factory->get_number_of_stored_characters())
                {
                    //std::cout << "Found using PARTIAL TEMP\n";
                    //std::cout << "Search hit during hybrid search\n";
                    return 1;
                }
                    
                cbc+=1;
                db>>=2;
                if (cbc >= 32){
                    cbc = 0;
                    dbi -= 1;
                    db = kmer_factory->blocks[dbi];
                }
                temp_array_block >>= 2;
                unchecked_temp_block_bits-=2;
                if (unchecked_temp_block_bits == 0)
                {
                    unchecked_temp_block_bits = 64;
                    current_temp_array_position -= 1;
                    temp_array_block = temp_array[current_temp_array_position];
                }
            }
        }
    }
    // If the characters have matched until now, the k-mer exists in the hash table
    return 1;
}

// Find where a k-mer is in the hash table
std::tuple<bool,uint64_t> PointerHashTable2::find_and_modify_using_prev(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = kmer_hash;
    // Do as long as probe sequence finds occupied hash table slots
    while (hash_table_array[kmer_slot].is_occupied())
    {
        // If previous k-mer exists, make a quick previous k-mer check
        if (previous_kmer_exists)
        {
            if (uint64_t(hash_table_array[kmer_slot].get_character()) == uint64_t(kmer_factory->get_newest_character()))
            {
                if (previous_kmer_slot == hash_table_array[kmer_slot].get_previous_kmer_slot())
                {
                    return std::make_tuple(true,kmer_slot);
                }
            }
        }
        // If quick check cannot be done or it is not successful, do a full check
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        // If k-mer in slot and the slot does not point to temp immediately
        if (check_result == 1)
        {
            return std::make_tuple(false,kmer_slot);
        }
        // If k-mer is not in the curent slot, probe to the next slot
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                resize();
                exit(1);
            }
        }
        // If k-mer is in the slot, byt the slot points to temp
        else if (check_result == 0)
        {
            // If we have a previous k-mer, migrate from temp to main
            if (previous_kmer_exists)
            {
                uint64_t temp_array_slot_in_question = hash_table_array[kmer_slot].get_previous_kmer_slot();
                // Release temp slot and recalculate
                temp_free_slots[temp_array_slot_in_question] = 1;
                temp_slots_in_use -= 1;
                for (int i = 0; i < kmer_factory->number_of_blocks; i++)
                {
                    temp_array[temp_array_slot_in_question*kmer_factory->number_of_blocks + i] = 0;
                }
                // Add main array count
                inserted_items+=1;
                hash_table_array[kmer_slot].set_prev_kmer_exists();
                hash_table_array[kmer_slot].set_previous_kmer_slot(previous_kmer_slot);   
                // Increase count
                //hash_table_array[kmer_slot].increase_count();
                // Set last character
                hash_table_array[kmer_slot].set_character(kmer_factory->get_newest_character());
                // Set occupied
                hash_table_array[kmer_slot].set_occupied();
                // Set in main array
                hash_table_array[kmer_slot].set_in_main_array();
                // Set complete/partial
                hash_table_array[kmer_slot].set_complete();
                inserted_complete_kmers += 1;
                // Set to not be at max count
                hash_table_array[kmer_slot].unset_at_max_count();

                return std::make_tuple(false,kmer_slot);
            }
            // If the new k-mer does not have a previous k-mer, do not migrate and instead increase count in temp
            else
            {
                //hash_table_array[kmer_slot].increase_count();
                return std::make_tuple(false,kmer_slot);
            }
        }
        else
        {
            std::cout << "ERROR : It should not be possible to end up in this situation\n";
            exit(1);
        }    
    }
    //std::cout << "SLOT " << kmer_slot <<  " WAS UNOCCUPIED\n";
    return std::make_tuple(false,size);
}

// Find where a k-mer is in the hash table
uint64_t PointerHashTable2::find(KMerFactory2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    while (hash_table_array[kmer_slot].is_occupied())
    {
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        if (check_result == 1)
        {
            return kmer_slot;
        } 
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                resize();
                exit(1);
            }
        }    
        else if (check_result == 0)
        {
            return kmer_slot;
        }
        else
        {
            std::cout << "ERROR : It should not be possible to end up in this situation\n";
            exit(1);
        }    
    }
    //std::cout << "SLOT " << kmer_slot <<  " WAS UNOCCUPIED\n";
    return size;
}

// Insert a new k-me to the hash table, we assume the k-mer is not already in the hash table
uint64_t PointerHashTable2::insert_new_kmer_in_main(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t prev_kmer_slot)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    
    while(hash_table_array[kmer_slot].is_occupied())
    {
        kmer_slot += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
        //kmer_slot += 1;
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == kmer_hash)
        {
            //std::cout << "Hash table is full, needs to be resized\n";
            resize();
            exit(1);
        }
    }

    // Next, add the k-mer info to the hash table

    // First, determine some characteristics of the k-mer
    bool complete = kmer_factory->get_number_of_stored_characters() >= kmer_len;
    bool prev_exists = kmer_factory->get_number_of_stored_characters() > 1;

    inserted_items+=1;

    // Set previous k-mer status (and previous k-mer slot) 
    if (prev_exists){
        hash_table_array[kmer_slot].set_prev_kmer_exists();
        hash_table_array[kmer_slot].set_previous_kmer_slot(prev_kmer_slot);
    } else {
        //std::cout << "DISASTROUS ERROR\n\n\n\n";
        hash_table_array[kmer_slot].unset_prev_kmer_exists();
    }
    // Set count to one
    hash_table_array[kmer_slot].set_count(1);
    // Set last character
    hash_table_array[kmer_slot].set_character(kmer_factory->get_newest_character());
    // Set occupied
    hash_table_array[kmer_slot].set_occupied();
    // Set in main array
    hash_table_array[kmer_slot].set_in_main_array();
    // Set complete/partial
    if (complete){
        hash_table_array[kmer_slot].set_complete();
        inserted_complete_kmers += 1;
    } else {
        hash_table_array[kmer_slot].unset_complete();
    }   
    // Set to not be at max count
    hash_table_array[kmer_slot].unset_at_max_count();
    //std::cout << "Inserted " << hash_table_array[kmer_slot].get_data() << "\n";
    //std::cout << "INSERTED IN MAIN USING SLOT " << kmer_slot << " IN MAIN AND PREV LINK TO SLOT " << prev_kmer_slot << "\n";
    return kmer_slot;
}

// Insert a new k-me to the hash table, we assume the k-mer is not already in the hash table
uint64_t PointerHashTable2::insert_new_kmer_in_temp(KMerFactory2BC* kmer_factory, uint64_t kmer_hash)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    
    while(hash_table_array[kmer_slot].is_occupied())
    {
        kmer_slot += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
        //kmer_slot += 1;
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == kmer_hash)
        {
            std::cout << "Hash table is full, needs to be resized\n";
            resize();
            exit(1);
        }
    }

    // Find the next free slot
    bool empty_temp_slot_found = false;
    for (int j = 0; j < max_temp_slots; j++)
    {
        if (temp_free_slots[j] == 1)
        {
            //std::cout << "Could be improved...\n";
            smallest_unused_temp_slot = j;
            empty_temp_slot_found = true;
            break;
        }
    }
    if (!empty_temp_slot_found)
    {
        std::cout << "TEMO SLOT RESIZING ERROR THAT SHOULD NOT HAPPEN\n";
        exit(1);
    }
    max_temp_slot_in_use = std::max(max_temp_slot_in_use, smallest_unused_temp_slot+1);
    temp_free_slots[smallest_unused_temp_slot] = 0;
    temp_slots_in_use += 1;
    
    hash_table_array[kmer_slot].set_occupied();
    hash_table_array[kmer_slot].unset_prev_kmer_exists();
    hash_table_array[kmer_slot].set_count(1);
    hash_table_array[kmer_slot].set_complete();
    hash_table_array[kmer_slot].unset_at_max_count();
    hash_table_array[kmer_slot].unset_in_main_array();
    hash_table_array[kmer_slot].set_previous_kmer_slot(smallest_unused_temp_slot);
    for (int i = 0; i < kmer_factory->number_of_blocks; i++)
    {
        temp_array[smallest_unused_temp_slot*kmer_factory->number_of_blocks + i] = kmer_factory->blocks[i];
        //std::cout << "WE PUT IN TEMP THE FOLLOWING K-MER: " << kmer_factory->blocks[i] << "\n";
    }

    if (smallest_unused_temp_slot == max_temp_slots-1)
    {
        std::cout << "Temp array resizing needed\n";
        std::cout << "Might be BUGGED.............................................\n";
        //exit(1);
        max_temp_slots *= 2;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
        //temp_free_slots = sdsl::bit_vector(max_temp_slots, 1);
        //temp_occupied_vector.resize(max_temp_slots);

    }
    //std::cout << "INSERTED IN TEMP USING SLOT " << kmer_slot << " IN MAIN AND SLOT " << used_temp_slots-1 << " IN TEMP\n";
    return kmer_slot;
}

uint64_t PointerHashTable2::find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_ip_step)
{
    bool using_prev;
    uint64_t kmer_slot;
    if (clean_ip_step)
    {
        kmer_slot = find(kmer_factory, kmer_hash);
    }
    else
    {
        std::tie(using_prev, kmer_slot) = find_and_modify_using_prev(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot);
    }
    if (kmer_slot < size)
    {
        if (kmer_factory->get_number_of_stored_characters() == kmer_len)
        {
            hash_table_array[kmer_slot].increase_count();
        }
        else
        {
            //std::cout << "===== SUS INCREASE TRY IN 'kmer_hash_table.cpp' function 'find_and_increment' =====\n";
        }
        
        if(previous_kmer_exists)
        {   
            if (hash_table_array[hash_table_array[kmer_slot].get_previous_kmer_slot()].get_count() < hash_table_array[previous_kmer_slot].get_count())
            {
                //std::cout << "PREVIOUS K-MER SWAPPED\n";
                hash_table_array[kmer_slot].set_previous_kmer_slot(previous_kmer_slot);
            }
        }
    }
        
    return kmer_slot;
}


uint64_t PointerHashTable2::find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_up_step)
{
    // Try to find and increment
    uint64_t kmer_slot = find_and_increment(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot, clean_up_step);
    // If that did not work, insert new
    if (kmer_slot == size)
    {
        // If previous k-mer exists, insert normally
        if (previous_kmer_exists || clean_up_step)
        {
            kmer_slot = insert_new_kmer_in_main(kmer_factory, kmer_hash, previous_kmer_slot);
            //std::cout << "Inserting new to main array\n";
        }
        // If previous k-mer does not exist, insert in temp
        else
        {
            kmer_slot = insert_new_kmer_in_temp(kmer_factory, kmer_hash);
            //std::cout << "Inserting new to temp array\n";
        }
    }
    else
    {
        //std::cout << "k-mer was found and its count increased\n";
    }
    return kmer_slot;
}

/*
uint64_t PointerHashTable2::clean_up_step_find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = clean_up_step_find_and_increment(kmer_factory, kmer_hash);
    if (kmer_slot == size)
    {
        kmer_slot = clean_up_step_insert_in_main(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot);
    }
    return kmer_slot;
}

uint64_t PointerHashTable2::clean_up_step_find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    while (hash_table_array[kmer_slot].is_occupied())
    {
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        if (check_result == 1)
        {
            hash_table_array[kmer_slot].increase_count();
            return kmer_slot;
        } 
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                resize();
                exit(1);
            }
        }    
        else if (check_result == 0)
        {
            std::cout << "ERROR!! THIS MUST NOT BE POSSIBLE (clean up step)\n";
            exit(1);
            hash_table_array[kmer_slot].increase_count();
            return kmer_slot;
        }
        else
        {
            std::cout << "ERROR!! It should not be possible to end up in this situation (clean up step)\n";
            exit(1);
        }    
    }
    return size;
}

uint64_t PointerHashTable2::clean_up_step_insert_in_main(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{

}
*/

void PointerHashTable2::clean_up_temp_array(KMerFactory2BC* kmer_factory, RollingHasher1* hasher)
{
    uint64_t brought_back_home = 0;
    for (int i = 0; i < size; i++)
    {
        if (!hash_table_array[i].is_occupied())
        {
            continue;
        }
        if (hash_table_array[i].is_in_main_array())
        {
            continue;
        }
        brought_back_home+=1;
        //std::cout << "K-MERS BROUGHT BACK HOME " << brought_back_home << "\n";
        //std::cout << "SLOT " << i << " IN MAIN WAS STORED IN TEMP, BRINGING IT BACK HOME\n";
        temp_slots_in_use -= 1;
        kmer_factory->reset();
        hasher->reset();

        uint64_t handled_characters = 0;
        uint64_t previous_kmer_slot = size;
        bool previous_kmer_exists = false;
        uint64_t slot_for_me = 0;
        uint64_t current_slot_in_temp = hash_table_array[i].get_previous_kmer_slot();

        while (handled_characters < kmer_len-1)
        {
            if (handled_characters > 0)
            {
                previous_kmer_exists = true;
            }
            kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
            hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            slot_for_me = hasher->get_current_hash();
            slot_for_me = find_and_increment_andifnot_insert(kmer_factory, slot_for_me, previous_kmer_exists, previous_kmer_slot, true);
            previous_kmer_slot = slot_for_me;
            handled_characters += 1;
        }
        slot_for_me = i;
        kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
        hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
        if (hasher->get_current_hash() != i)
        {
            //std::cout << "SUS ERROR !!\n";
        }

        hash_table_array[slot_for_me].set_prev_kmer_exists();
        hash_table_array[slot_for_me].set_previous_kmer_slot(previous_kmer_slot);
        hash_table_array[slot_for_me].set_character(kmer_factory->get_newest_character());
        hash_table_array[slot_for_me].set_in_main_array();
        hash_table_array[slot_for_me].set_complete();
        inserted_complete_kmers += 1;
        inserted_items+=1;
    }
}


void PointerHashTable2::clean_up_temp_array_OLD(KMerFactory2BC* kmer_factory, RollingHasher1* hasher)
{
    uint64_t brought_back_home = 0;
    for (int i = 0; i <= size; i++)
    {
        if (!hash_table_array[i].is_occupied())
        {
            continue;
        }
        if (hash_table_array[i].is_in_main_array())
        {
            continue;
        }
        brought_back_home+=1;
        //std::cout << "K-MERS BROUGHT BACK HOME " << brought_back_home << "\n";
        //std::cout << "SLOT " << i << " IN MAIN WAS STORED IN TEMP, BRINGING IT BACK HOME\n";
        temp_slots_in_use -= 1;
        kmer_factory->reset();
        hasher->reset();

        uint64_t handled_characters = 0;
        uint64_t previous_kmer_slot = size;
        uint64_t slot_for_me = 0;
        uint64_t current_slot_in_temp = hash_table_array[i].get_previous_kmer_slot();

        while (handled_characters < kmer_len-1)
        {
            kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
            hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            slot_for_me = hasher->get_current_hash();
            while(hash_table_array[slot_for_me].is_occupied())
            {
                slot_for_me += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
                if (slot_for_me >= size)
                    slot_for_me = slot_for_me % size;
                if (slot_for_me == hasher->get_current_hash())
                {
                    std::cout << "Hash table is full, needs to be resized\n";
                    resize();
                    exit(1);
                }
            }
            insert_new_kmer_in_main(kmer_factory, slot_for_me, previous_kmer_slot);
            previous_kmer_slot = slot_for_me;
            handled_characters += 1;
            //std::cout << "HANDLED CHARACTERS " << handled_characters << "\n";
        }
        slot_for_me = i;
        kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
        hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
        if (hasher->get_current_hash() != i)
        {
            //std::cout << "SUS ERROR !!\n";
        }
        hash_table_array[slot_for_me].set_prev_kmer_exists();
        hash_table_array[slot_for_me].set_previous_kmer_slot(previous_kmer_slot);
        hash_table_array[slot_for_me].set_character(kmer_factory->get_newest_character());
        hash_table_array[slot_for_me].set_in_main_array();
        hash_table_array[slot_for_me].set_complete();
        inserted_complete_kmers += 1;
        inserted_items+=1;
    }
}


char PointerHashTable2::get_oldest_char_in_temp_slot_and_shift_left(uint64_t slot, KMerFactory2BC* kmer_factory)
{
    uint64_t next_char_chunk = ((temp_array[slot*kmer_factory->number_of_blocks]) >> (kmer_factory->bits_in_last_block-2));
    char next_char = twobitstringfunctions::int2char(next_char_chunk);
    for (int i = 0; i < kmer_factory->number_of_blocks-1; i++)
    {
        temp_array[slot*kmer_factory->number_of_blocks + i] = temp_array[slot*kmer_factory->number_of_blocks + i] << 2;
        temp_array[slot*kmer_factory->number_of_blocks + i] = (temp_array[slot*kmer_factory->number_of_blocks + i] | (temp_array[slot*kmer_factory->number_of_blocks + i + 1] >> 64 - 2));
    }
    temp_array[slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1] = (temp_array[slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1] << 2);
    temp_array[slot*kmer_factory->number_of_blocks] = (temp_array[slot*kmer_factory->number_of_blocks] & (kmer_factory->used_left_block_mask));
    return next_char;
}

uint64_t PointerHashTable2::get_number_of_max_temp_slots()
{
    return max_temp_slots;
}

uint64_t PointerHashTable2::get_number_of_touched_temp_slots()
{
    return touched_temp_slots;
}

uint64_t PointerHashTable2::get_number_of_temp_slots_in_use()
{
    return temp_slots_in_use;
}

uint64_t PointerHashTable2::get_number_of_max_temp_slots_in_use()
{
    return max_temp_slot_in_use;
}

uint64_t PointerHashTable2::get_smallest_unused_temp_slot()
{
    for (int j = 0; j < max_temp_slots; j++)
    {
        if (temp_free_slots[j] == 1)
        {
            smallest_unused_temp_slot = j;
            break;
        }
    }
    return smallest_unused_temp_slot;
}

void PointerHashTable2::write_kmers_on_disk(KMerFactory2BC* kmer_factory, uint64_t min_abundance, bool only_canonical, std::string output_path)
{
    //std::string output_file_path = output_directory + "/" + std::to_string(kmer_len) + "-mers_" + std::to_string(min_abundance) + "-abundance_2-mode.txt";
	std::ofstream output_file(output_path);

    uint64_t check_position = 0;
    uint64_t chain_position = 0;
    std::vector<uint64_t> current_chain_abundances;
    //std::vector<bool> current_chain_completeness;
    std::string current_chain_string = "";
    int current_chain_streak = 0;
    int round_one_kmers = 0;
    int round_two_kmers = 0;

    sdsl::bit_vector is_referenced(size,0);
    while(check_position < size)
    {
        if (hash_table_array[check_position].is_occupied())
        {
            is_referenced[hash_table_array[check_position].get_previous_kmer_slot()] = 1;
        }
        check_position+=1;
    }

    for (int round = 0; round < 2; round++)
    {
        //if (round == 0)
        //{
        //    continue;
        //}          
        check_position = 0;
        chain_position = 0;
        current_chain_abundances.clear();
        current_chain_string = "";
        current_chain_streak = 0;

        while (check_position < size)
        {
            if (round == 0)
            {
                if (is_referenced[check_position] == 1)
                {
                    check_position+=1;
                    continue;
                }
            }
            if ((hash_table_array[check_position].is_complete()) && (!hash_table_array[check_position].is_written_in_output()) && (hash_table_array[check_position].is_occupied()) && (hash_table_array[check_position].is_complete()) && (hash_table_array[check_position].get_count() >= min_abundance))
            {
                chain_position = check_position;
                // Chain together all k-mers that have not been written into output
                while (!hash_table_array[chain_position].is_written_in_output() && (hash_table_array[chain_position].is_complete()))
                {
                    current_chain_streak+=1;
                    current_chain_abundances.push_back(hash_table_array[chain_position].get_count());
                    //current_chain_completeness.push_back(hash_table_array[chain_position].is_complete());
                    current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character())) + current_chain_string;
                    hash_table_array[chain_position].set_is_written_in_output();
                    if(hash_table_array[chain_position].prev_kmer_exists())
                    {
                        chain_position = hash_table_array[chain_position].get_previous_kmer_slot();
                    }
                    else
                    {
                        break;
                    }
                }
                // Complete the last k-mer in the chain
                for (int i = 0; i < kmer_len-1; i++)
                {
                    current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character())) + current_chain_string;
                    if (hash_table_array[chain_position].prev_kmer_exists())
                    {
                        chain_position = hash_table_array[chain_position].get_previous_kmer_slot();
                    }
                    else
                    {
                        break;
                    }
                }
                // Write the chained k-mer into output
                for (int j = 0; j < current_chain_streak; j++)
                {
                    if (j+kmer_len > current_chain_string.length())
                    {
                        break;
                    }
                    if (round == 0)
                        round_one_kmers+=1;
                    else
                        round_two_kmers+=1;

                    std::string current_kmer = current_chain_string.substr(j, kmer_len);
                    int current_kmer_is_canonical = purestringfunctions::is_canonical(current_kmer);

                    if (current_kmer_is_canonical == 1)
                    {
                        if ((current_chain_abundances[current_chain_streak-j-1] >= min_abundance))
                        {
                            output_file << current_chain_string.substr(j, kmer_len) << " " << std::to_string(current_chain_abundances[current_chain_streak-j-1]) << "\n";
                        }
                    }
                    else if (current_kmer_is_canonical == 0)
                    {
                        if ((current_chain_abundances[current_chain_streak-j-1]/2 >= min_abundance))
                        {
                            output_file << current_chain_string.substr(j, kmer_len) << " " << std::to_string(current_chain_abundances[current_chain_streak-j-1]/2) << "\n";
                        }
                    }
                }
                // Reset 
                //current_chain_completeness.clear();
                current_chain_abundances.clear();
                current_chain_string = "";
                current_chain_streak = 0;
            }
            check_position += 1;
        }
    }
    std::cout << "k-mer starts in round one: " << round_one_kmers << "\n";
    std::cout << "k-mer starts in round two: " << round_two_kmers << "\n";

    output_file.close(); // close file
	output_file.clear(); // clear flags
} 

uint64_t PointerHashTable2::get_size()
{
    return size;
}

void PointerHashTable2::copy_content_from_another_hash_table(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table)
{
    if (max_temp_slots != old_hash_table->max_temp_slots)
    {
        max_temp_slots = old_hash_table->max_temp_slots;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
    }
    // Copy content from temp array
    for (int i = 0; i < old_hash_table->temp_array.size(); i++)
    {
        this->temp_array[i] = old_hash_table->temp_array[i];
    }
    // Also remember to copy the bit vector here
    for (int i = 0; i < old_hash_table->temp_free_slots.size(); i++)
    {
        this->temp_free_slots[i] = old_hash_table->temp_free_slots[i];
    }
    inserted_items = old_hash_table->inserted_items;
    inserted_complete_kmers = old_hash_table->inserted_complete_kmers;
    touched_temp_slots = old_hash_table->touched_temp_slots;
    temp_slots_in_use = old_hash_table->temp_slots_in_use;
    smallest_unused_temp_slot = old_hash_table->smallest_unused_temp_slot;

    

    // Start checking slots in the old hash table
    uint64_t check_position = 0;
    uint64_t chain_position;
    uint64_t chain_length;
    std::vector<uint64_t> old_positions_chain;
    //std::vector<uint64_t> where_did_i_go(old_hash_table->size, 0);
    int percentage_checked = 0;
    uint64_t migrated_count = 0;
    
    sdsl::bit_vector is_referenced(old_hash_table->get_size(),0);
    while(check_position < old_hash_table->get_size())
    {
        if (old_hash_table->hash_table_array[check_position].is_occupied())
        {
            is_referenced[old_hash_table->hash_table_array[check_position].get_previous_kmer_slot()] = 1;
        }
        check_position+=1;
    }

    int round_one_kmers = 0;
    int round_two_kmers = 0;

    // Flag 1 used in loop 
    // FLag 2 used to check if already transferred
    int percentage = 1;

    for (int round = 0; round < 2; round++)
    {
        check_position = 0;

        while (check_position < old_hash_table->get_size())
        {
            if (round == 0)
            {
                if (is_referenced[check_position] == 1)
                {
                    check_position+=1;
                    continue;
                }
            }

            //std::cout << "CHECK POSITION " << check_position << "\n";
            if (((100.0*migrated_count) / old_hash_table->get_size()) >=percentage)
            {
                std::cout << "Resizing... " << percentage << "%" << " done\n";
                percentage+=1;
            }
            //if (check_position % 100000 == 0)
            //    std::cout << "Resizing... " << check_position << " / " << old_hash_table->size << "\n";
            // If the k-mer at the check position is already migrated, continue to next position
            if (old_hash_table->hash_table_array[check_position].is_flagged_2())
            {
                //std::cout << "## k-mer already migrated at position " << check_position << "\n";
                check_position += 1;
                continue;
            }
            if (!old_hash_table->hash_table_array[check_position].is_occupied())
            {
                //std::cout << "?? no k-mer at position " << check_position << "\n";
                old_hash_table->hash_table_array[check_position].set_is_flagged_2();
                check_position += 1;
                migrated_count+=1;
                continue;
            }
            // Reset the previous chain
            old_positions_chain.clear();
            kmer_factory->reset();
            new_hasher->reset();

            // If the k-mer has not yet been migrated, create the chain and migrate all
            chain_position = check_position;
            chain_length = 0;
            // Chain all k-mers that have not been migrated AND that are not laready in the loop
            bool break_out_because_no_previous = false;
            while((!old_hash_table->hash_table_array[chain_position].is_flagged_1()) && (!old_hash_table->hash_table_array[chain_position].is_flagged_2()))
            {
                //std::cout << "CHAIN POSITION IS " << chain_position << "\n";
                //std::cout << "- its previous k-mer is at position " << old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot() << "\n";
                // Add the chain position to the chain position list
                old_positions_chain.push_back(chain_position);
                chain_length += 1;
                // Mark the k-mer so that we do not include it twice in the loop
                old_hash_table->hash_table_array[chain_position].set_is_flagged_1();
                // If previous k-mer exists, we check it next
                if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                {
                    if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                    {
                        chain_position = hash_table_array[old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot()].get_previous_kmer_slot();
                        std::cout << "SUS BEHAVIOR IF WE GET HERE\n";
                    }
                    else
                    {
                        chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                    }
                }
                // Otherwise we leave the loop
                else
                {
                    //std::cout << "~~ this position does not have previous k-mer\n";
                    break_out_because_no_previous = true;
                    break;
                }
            }
            // Unset loop flags
            for (uint64_t op : old_positions_chain)
            {
                old_hash_table->hash_table_array[op].unset_is_flagged_1();
            }
            // Then add k-1 previous k-mers so we can calculate hash values for all k-mers
            if (!break_out_because_no_previous)
            {
                //std::cout << "ADDING EXTRA K-MERS OT COMPLETE THE LAST ONE\n";
                for (int j = 0; j < kmer_len-1; j++)
                {
                    if (old_hash_table->hash_table_array[chain_position].is_occupied())
                    {
                        old_positions_chain.push_back(chain_position);
                        chain_length += 1;
                        //std::cout << "ALSO ADDED CHAIN POSITION " << chain_position << "\n";
                    }
                    else
                    {
                        //std::cout << "COULD NOT ADD CHAIN POSITION " << chain_position << "\n";
                        break;
                    }
                    // If the previous k-mer exists, we add it to the list
                    if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                    {
                        if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                        {
                            chain_position = hash_table_array[old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot()].get_previous_kmer_slot();
                        }
                        else
                        {
                            chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                        }
                        //std::cout << "this one has a previous k-mer at " << chain_position << "\n";
                    }
                    // Otherwise break the loop early
                    else
                    {
                        //std::cout << "this one did not have previous k-mer\n";
                        //std::cout << "NO PREVIOUS K-MER EXISTS IN CAHIN POSITION " << chain_position << "\n";
                        if (old_hash_table->hash_table_array[chain_position].is_in_main_array())
                        {
                            //std::cout << "BUT IS IN MAIN ARRAY??????\n";
                        }
                        break;
                    }
                }
            }
            

            // Now, start pushing characters from the chain (in reverse order) into the k-mer factory and calculate their new hash values
            
            // First, check if the last k-mer is in the temp array
            //std::cout << "FINDING FULL CHAIN for chain where initial length is: " << chain_length << "\n";
            uint64_t at_chain_position = chain_length-1;
            //std::cout << "THERE ARE THIS MANY ITEMS IN THE CHAIN VECTOR " << old_positions_chain.size() << "\n";
            uint64_t last_position_in_chain = old_positions_chain[at_chain_position];
            //std::cout << "The last position in chain is: " << last_position_in_chain << "\n";
            // If the last k-mer is in temp, load it fully into k-mer factory
            uint64_t old_temp_slot;
            if (!old_hash_table->hash_table_array[last_position_in_chain].is_in_main_array())
            {
                //std::cout << "was not in main\n";
                // Get the slot in temp
                if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                {
                    old_temp_slot = hash_table_array[old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot()].get_previous_kmer_slot();
                }
                else
                {
                    old_temp_slot = old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot();
                }
               
                for (int k = 0; k < kmer_factory->number_of_blocks; k++)
                {
                    //std::cout << "POSSIBLE ERROR SPOT STARTS\n";
                    uint64_t block_to_be_pushed = old_hash_table->temp_array[old_temp_slot*kmer_blocks+k];
                    //std::cout << "GOT OVER IT\n";
                    // If we are pushing the first block do this (possibly not full block)
                    if (k == 0)
                    {
                        //std::cout << "...\n";
                        for (int m = 0; m < kmer_factory->bits_in_last_block; m+=2)
                        {
                            //std::cout << "pushed char: " << (uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2))) << "\n"; 
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                            
                        }
                        //std::cout << "...\n";
                    }
                    // Otherwise we are pushing a full block
                    else 
                    {
                        for (int m = 0; m < 64; m+=2)
                        {
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (64-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                        }
                    }
                }
            } // Now the kmer_factory holds the "first" k-mer completely
            // Alternatively, if the last k-mer is in the main array do this
            else
            {
                //std::cout << "was in main\n";
                at_chain_position = chain_length-1;
                //std::cout << "...\n";
                for (int k = 0; k < kmer_len; k++)
                {
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    at_chain_position-=1;
                }
                //std::cout << "...\n";
                at_chain_position+=1;
            }

            //std::cout << "FULL CHAIN FOUND\n";
            //std::cout << "START ADDING FROM OLD POSITION " << old_positions_chain[at_chain_position] << "\n";

            // Now k-mer is loaded with a k-mer whose position in the old hash table is in chain_positions[at_chain_position]

            // Start checking the k-mers and add them to the new hash table it flag2 is not set
            uint64_t at_chain_position_backup = at_chain_position;
            while(at_chain_position >= 0)
            {
                // Get the position in the old array
                uint64_t old_position = old_positions_chain[at_chain_position];
                //std::cout << "- MIGRATING k-MER FROM OLD POSITION " << old_position << "\n";
                
                // If not flagged2, migrate to new hash table
                if (!old_hash_table->hash_table_array[old_position].is_flagged_2())
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << old_hash_table->hash_table_array[old_position].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- actually transfering position " << old_position << "\n";
                    // Initial new position
                    uint64_t new_position = new_hasher->get_current_hash();
                    uint64_t probing_backup = new_position;
                    // Probe for empty new position
                    while(hash_table_array[new_position].is_occupied())
                    {
                        new_position += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
                        if (new_position >= size)
                            new_position = new_position % size;
                        if (new_position == probing_backup){
                            std::cout << "RESIZING FAILURE, VERY LETHAL\n";
                            exit(1);
                        }
                    }

                    // Insert into the found empty position
                    hash_table_array[new_position].set_data(old_hash_table->hash_table_array[old_position].get_data());
                    hash_table_array[new_position].set_count(old_hash_table->hash_table_array[old_position].get_count());
                    //hash_table_array[new_position].set_previous_kmer_slot(old_position);
                    hash_table_array[new_position].set_previous_kmer_slot(old_hash_table->hash_table_array[old_position].get_previous_kmer_slot());
                    hash_table_array[new_position].unset_is_flagged_1();
                    hash_table_array[new_position].unset_is_flagged_2();
                    hash_table_array[new_position].unset_is_written_in_output();
                    
                    // Remember where the old position was migrated to
                    //where_did_i_go[old_position] = new_position;
                    old_hash_table->hash_table_array[old_position].set_previous_kmer_slot(new_position);
                    // Set flag
                    old_hash_table->hash_table_array[old_position].set_is_flagged_2();
                    migrated_count+=1;

                    if (round == 0)
                        round_one_kmers+=1;
                    else
                        round_two_kmers+=1;

                    // Count new slot
                    // Add
                    // Flags
                }
                else
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << hash_table_array[old_hash_table->hash_table_array[old_position].get_previous_kmer_slot()].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- not transfering position " << old_position <<  " because already transferred\n";
                }

                at_chain_position-=1;
                // All is done
                if (at_chain_position == -1)
                {
                    break;
                }
                // Or continue with the next
                else
                {
                    //std::cout << "...\n";
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    //std::cout << "...\n";
                }
            }
            //old_hash_table->hash_table_array[check_position].set_is_flagged_2();
            check_position+=1;
        }
    }

    // SET THE CORRECT PREVIOUS K-MERS IN THE NEW HASH TABLE
    uint64_t check_position_old;
    uint64_t check_position_old_previous;
    uint64_t check_position_new_previous;
    check_position = 0;
    while (check_position < size)
    {
        // If not occupied, skip
        if (!hash_table_array[check_position].is_occupied())
        {
            check_position += 1;
            continue;
        }
        // If in main array, make the swap
        if (hash_table_array[check_position].is_in_main_array())
        {
            hash_table_array[check_position].set_previous_kmer_slot(old_hash_table->hash_table_array[hash_table_array[check_position].get_previous_kmer_slot()].get_previous_kmer_slot());
        }
        // Otherwise, swap is not necessary
        /*
        hash_table_array[check_position].set_previous_kmer_slot();
        check_position_old = hash_table_array[check_position].get_previous_kmer_slot();
        // If old k-mer was in main
        if (old_hash_table->hash_table_array[check_position_old].is_in_main_array())
        {
            check_position_old_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
            check_position_new_previous = where_did_i_go[check_position_old_previous];
        }
        // And if it was not
        else
        {
            check_position_new_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
        }
        hash_table_array[check_position].set_previous_kmer_slot(check_position_new_previous);
        }*/
        
        check_position += 1;
    }
    std::cout << "Resizing finished\n";
    std::cout << "k-mers migrated in round one: " << round_one_kmers << "\n";
    std::cout << "k-mers migrated in round two: " << round_two_kmers << "\n";
}

void PointerHashTable2::copy_content_from_another_hash_table_BIGMEM(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table)
{
    if (max_temp_slots != old_hash_table->max_temp_slots)
    {
        max_temp_slots = old_hash_table->max_temp_slots;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
    }
    // Copy content from temp array
    for (int i = 0; i < old_hash_table->temp_array.size(); i++)
    {
        this->temp_array[i] = old_hash_table->temp_array[i];
    }
    // Also remember to copy the bit vector here
    for (int i = 0; i < old_hash_table->temp_free_slots.size(); i++)
    {
        this->temp_free_slots[i] = old_hash_table->temp_free_slots[i];
    }
    inserted_items = old_hash_table->inserted_items;
    inserted_complete_kmers = old_hash_table->inserted_complete_kmers;
    touched_temp_slots = old_hash_table->touched_temp_slots;
    temp_slots_in_use = old_hash_table->temp_slots_in_use;
    smallest_unused_temp_slot = old_hash_table->smallest_unused_temp_slot;

    

    // Start checking slots in the old hash table
    uint64_t check_position = 0;
    uint64_t chain_position;
    uint64_t chain_length;
    std::vector<uint64_t> old_positions_chain;
    std::vector<uint64_t> where_did_i_go(old_hash_table->size, 0);
    int percentage_checked = 0;
    uint64_t migrated_count = 0;
    
    sdsl::bit_vector is_referenced(old_hash_table->get_size(),0);
    while(check_position < old_hash_table->get_size())
    {
        if (old_hash_table->hash_table_array[check_position].is_occupied())
        {
            is_referenced[old_hash_table->hash_table_array[check_position].get_previous_kmer_slot()] = 1;
        }
        check_position+=1;
    }

    int round_one_kmers = 0;
    int round_two_kmers = 0;

    // Flag 1 used in loop 
    // FLag 2 used to check if already transferred
    int percentage = 1;

    for (int round = 0; round < 2; round++)
    {
        check_position = 0;

        while (check_position < old_hash_table->get_size())
        {
            if (round == 0)
            {
                if (is_referenced[check_position] == 1)
                {
                    check_position+=1;
                    continue;
                }
            }

            //std::cout << "CHECK POSITION " << check_position << "\n";
            if (((100.0*migrated_count) / old_hash_table->get_size()) >=percentage)
            {
                std::cout << "Resizing... " << percentage << "%" << " done\n";
                percentage+=1;
            }
            //if (check_position % 100000 == 0)
            //    std::cout << "Resizing... " << check_position << " / " << old_hash_table->size << "\n";
            // If the k-mer at the check position is already migrated, continue to next position
            if (old_hash_table->hash_table_array[check_position].is_flagged_2())
            {
                //std::cout << "## k-mer already migrated at position " << check_position << "\n";
                check_position += 1;
                continue;
            }
            if (!old_hash_table->hash_table_array[check_position].is_occupied())
            {
                //std::cout << "?? no k-mer at position " << check_position << "\n";
                old_hash_table->hash_table_array[check_position].set_is_flagged_2();
                check_position += 1;
                migrated_count+=1;
                continue;
            }
            // Reset the previous chain
            old_positions_chain.clear();
            kmer_factory->reset();
            new_hasher->reset();

            // If the k-mer has not yet been migrated, create the chain and migrate all
            chain_position = check_position;
            chain_length = 0;
            // Chain all k-mers that have not been migrated AND that are not laready in the loop
            bool break_out_because_no_previous = false;
            while((!old_hash_table->hash_table_array[chain_position].is_flagged_1()) && (!old_hash_table->hash_table_array[chain_position].is_flagged_2()))
            {
                //std::cout << "CHAIN POSITION IS " << chain_position << "\n";
                //std::cout << "- its previous k-mer is at position " << old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot() << "\n";
                // Add the chain position to the chain position list
                old_positions_chain.push_back(chain_position);
                chain_length += 1;
                // Mark the k-mer so that we do not include it twice in the loop
                old_hash_table->hash_table_array[chain_position].set_is_flagged_1();
                // If previous k-mer exists, we check it next
                if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                {
                    chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                }
                // Otherwise we leave the loop
                else
                {
                    //std::cout << "~~ this position does not have previous k-mer\n";
                    break_out_because_no_previous = true;
                    break;
                }
            }
            // Unset loop flags
            for (uint64_t op : old_positions_chain)
            {
                old_hash_table->hash_table_array[op].unset_is_flagged_1();
            }
            // Then add k-1 previous k-mers so we can calculate hash values for all k-mers
            if (!break_out_because_no_previous)
            {
                //std::cout << "ADDING EXTRA K-MERS OT COMPLETE THE LAST ONE\n";
                for (int j = 0; j < kmer_len-1; j++)
                {
                    if (old_hash_table->hash_table_array[chain_position].is_occupied())
                    {
                        old_positions_chain.push_back(chain_position);
                        chain_length += 1;
                        //std::cout << "ALSO ADDED CHAIN POSITION " << chain_position << "\n";
                    }
                    else
                    {
                        //std::cout << "COULD NOT ADD CHAIN POSITION " << chain_position << "\n";
                        break;
                    }
                    // If the previous k-mer exists, we add it to the list
                    if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                    {
                        chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                        //std::cout << "this one has a previous k-mer at " << chain_position << "\n";
                    }
                    // Otherwise break the loop early
                    else
                    {
                        //std::cout << "this one did not have previous k-mer\n";
                        //std::cout << "NO PREVIOUS K-MER EXISTS IN CAHIN POSITION " << chain_position << "\n";
                        if (old_hash_table->hash_table_array[chain_position].is_in_main_array())
                        {
                            //std::cout << "BUT IS IN MAIN ARRAY??????\n";
                        }
                        break;
                    }
                }
            }
            

            // Now, start pushing characters from the chain (in reverse order) into the k-mer factory and calculate their new hash values
            
            // First, check if the last k-mer is in the temp array
            //std::cout << "FINDING FULL CHAIN for chain where initial length is: " << chain_length << "\n";
            uint64_t at_chain_position = chain_length-1;
            //std::cout << "THERE ARE THIS MANY ITEMS IN THE CHAIN VECTOR " << old_positions_chain.size() << "\n";
            uint64_t last_position_in_chain = old_positions_chain[at_chain_position];
            //std::cout << "The last position in chain is: " << last_position_in_chain << "\n";
            // If the last k-mer is in temp, load it fully into k-mer factory
            if (!old_hash_table->hash_table_array[last_position_in_chain].is_in_main_array())
            {
                //std::cout << "was not in main\n";
                // Get the slot in temp
                uint64_t old_temp_slot = old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot();
                for (int k = 0; k < kmer_factory->number_of_blocks; k++)
                {
                    //std::cout << "POSSIBLE ERROR SPOT STARTS\n";
                    uint64_t block_to_be_pushed = old_hash_table->temp_array[old_temp_slot*kmer_blocks+k];
                    //std::cout << "GOT OVER IT\n";
                    // If we are pushing the first block do this (possibly not full block)
                    if (k == 0)
                    {
                        //std::cout << "...\n";
                        for (int m = 0; m < kmer_factory->bits_in_last_block; m+=2)
                        {
                            //std::cout << "pushed char: " << (uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2))) << "\n"; 
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                            
                        }
                        //std::cout << "...\n";
                    }
                    // Otherwise we are pushing a full block
                    else 
                    {
                        for (int m = 0; m < 64; m+=2)
                        {
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (64-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                        }
                    }
                }
            } // Now the kmer_factory holds the "first" k-mer completely
            // Alternatively, if the last k-mer is in the main array do this
            else
            {
                //std::cout << "was in main\n";
                at_chain_position = chain_length-1;
                //std::cout << "...\n";
                for (int k = 0; k < kmer_len; k++)
                {
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    at_chain_position-=1;
                }
                //std::cout << "...\n";
                at_chain_position+=1;
            }

            //std::cout << "FULL CHAIN FOUND\n";
            //std::cout << "START ADDING FROM OLD POSITION " << old_positions_chain[at_chain_position] << "\n";

            // Now k-mer is loaded with a k-mer whose position in the old hash table is in chain_positions[at_chain_position]

            // Start checking the k-mers and add them to the new hash table it flag2 is not set
            uint64_t at_chain_position_backup = at_chain_position;
            while(at_chain_position >= 0)
            {
                // Get the position in the old array
                uint64_t old_position = old_positions_chain[at_chain_position];
                //std::cout << "- MIGRATING k-MER FROM OLD POSITION " << old_position << "\n";
                
                // If not flagged2, migrate to new hash table
                if (!old_hash_table->hash_table_array[old_position].is_flagged_2())
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << old_hash_table->hash_table_array[old_position].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- actually transfering position " << old_position << "\n";
                    // Initial new position
                    uint64_t new_position = new_hasher->get_current_hash();
                    uint64_t probing_backup = new_position;
                    // Probe for empty new position
                    while(hash_table_array[new_position].is_occupied())
                    {
                        new_position += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
                        if (new_position >= size)
                            new_position = new_position % size;
                        if (new_position == probing_backup){
                            std::cout << "RESIZING FAILURE, VERY LETHAL\n";
                            exit(1);
                        }
                    }

                    // Insert into the found empty position
                    hash_table_array[new_position].set_data(old_hash_table->hash_table_array[old_position].get_data());
                    hash_table_array[new_position].set_count(old_hash_table->hash_table_array[old_position].get_count());
                    hash_table_array[new_position].set_previous_kmer_slot(old_position);
                    hash_table_array[new_position].unset_is_flagged_1();
                    hash_table_array[new_position].unset_is_flagged_2();
                    hash_table_array[new_position].unset_is_written_in_output();
                    
                    // Remember where the old position was migrated to
                    where_did_i_go[old_position] = new_position;
                    // Set flag
                    old_hash_table->hash_table_array[old_position].set_is_flagged_2();
                    migrated_count+=1;

                    if (round == 0)
                        round_one_kmers+=1;
                    else
                        round_two_kmers+=1;

                    // Count new slot
                    // Add
                    // Flags
                }
                else
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << hash_table_array[old_hash_table->hash_table_array[old_position].get_previous_kmer_slot()].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- not transfering position " << old_position <<  " because already transferred\n";
                }

                at_chain_position-=1;
                // All is done
                if (at_chain_position == -1)
                {
                    break;
                }
                // Or continue with the next
                else
                {
                    //std::cout << "...\n";
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    //std::cout << "...\n";
                }
            }
            //old_hash_table->hash_table_array[check_position].set_is_flagged_2();
            check_position+=1;
        }
    }

    // SET THE CORRECT PREVIOUS K-MERS IN THE NEW HASH TABLE
    uint64_t check_position_old;
    uint64_t check_position_old_previous;
    uint64_t check_position_new_previous;
    check_position = 0;
    while (check_position < size)
    {
        if (hash_table_array[check_position].is_occupied())
        {
            check_position_old = hash_table_array[check_position].get_previous_kmer_slot();
            // If old k-mer was in main
            if (old_hash_table->hash_table_array[check_position_old].is_in_main_array())
            {
                check_position_old_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
                check_position_new_previous = where_did_i_go[check_position_old_previous];
            }
            // And if it was not
            else
            {
                check_position_new_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
            }
            
            hash_table_array[check_position].set_previous_kmer_slot(check_position_new_previous);
        }
        check_position += 1;
    }
    std::cout << "Resizing finished\n";
    std::cout << "k-mers migrated in round one: " << round_one_kmers << "\n";
    std::cout << "k-mers migrated in round two: " << round_two_kmers << "\n";
}


void PointerHashTable2::copy_content_from_another_hash_table_OLD_AND_WORKING(KMerFactory2BC* kmer_factory, RollingHasher1* new_hasher, PointerHashTable2* old_hash_table)
{
    if (max_temp_slots != old_hash_table->max_temp_slots)
    {
        max_temp_slots = old_hash_table->max_temp_slots;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
    }
    // Copy content from temp array
    for (int i = 0; i < old_hash_table->temp_array.size(); i++)
    {
        this->temp_array[i] = old_hash_table->temp_array[i];
    }
    // Also remember to copy the bit vector here
    for (int i = 0; i < old_hash_table->temp_free_slots.size(); i++)
    {
        this->temp_free_slots[i] = old_hash_table->temp_free_slots[i];
    }
    inserted_items = old_hash_table->inserted_items;
    inserted_complete_kmers = old_hash_table->inserted_complete_kmers;
    touched_temp_slots = old_hash_table->touched_temp_slots;
    temp_slots_in_use = old_hash_table->temp_slots_in_use;
    smallest_unused_temp_slot = old_hash_table->smallest_unused_temp_slot;

    

    // Start checking slots in the old hash table
    uint64_t check_position = 0;
    uint64_t chain_position;
    uint64_t chain_length;
    std::vector<uint64_t> old_positions_chain;
    std::vector<uint64_t> where_did_i_go(old_hash_table->size, 0);
    int percentage_checked = 0;
    uint64_t migrated_count = 0;
    
    // Flag 1 used in loop 
    // FLag 2 used to check if already transferred
    int percentage = 1;
    while (check_position < old_hash_table->get_size())
    {
        //std::cout << "CHECK POSITION " << check_position << "\n";
        if (((100.0*migrated_count) / old_hash_table->get_size()) >=percentage)
        {
            std::cout << "Resizing... " << percentage << "%" << " done\n";
            percentage+=1;
        }
        //if (check_position % 100000 == 0)
        //    std::cout << "Resizing... " << check_position << " / " << old_hash_table->size << "\n";
        // If the k-mer at the check position is already migrated, continue to next position
        if (old_hash_table->hash_table_array[check_position].is_flagged_2())
        {
            //std::cout << "## k-mer already migrated at position " << check_position << "\n";
            check_position += 1;
            continue;
        }
        if (!old_hash_table->hash_table_array[check_position].is_occupied())
        {
            //std::cout << "?? no k-mer at position " << check_position << "\n";
            old_hash_table->hash_table_array[check_position].set_is_flagged_2();
            check_position += 1;
            migrated_count+=1;
            continue;
        }
        // Reset the previous chain
        old_positions_chain.clear();
        kmer_factory->reset();
        new_hasher->reset();

        // If the k-mer has not yet been migrated, create the chain and migrate all
        chain_position = check_position;
        chain_length = 0;
        // Chain all k-mers that have not been migrated AND that are not laready in the loop
        bool break_out_because_no_previous = false;
        while((!old_hash_table->hash_table_array[chain_position].is_flagged_1()) && (!old_hash_table->hash_table_array[chain_position].is_flagged_2()))
        {
            //std::cout << "CHAIN POSITION IS " << chain_position << "\n";
            //std::cout << "- its previous k-mer is at position " << old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot() << "\n";
            // Add the chain position to the chain position list
            old_positions_chain.push_back(chain_position);
            chain_length += 1;
            // Mark the k-mer so that we do not include it twice in the loop
            old_hash_table->hash_table_array[chain_position].set_is_flagged_1();
            // If previous k-mer exists, we check it next
            if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
            {
                chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
            }
            // Otherwise we leave the loop
            else
            {
                //std::cout << "~~ this position does not have previous k-mer\n";
                break_out_because_no_previous = true;
                break;
            }
        }
        // Unset loop flags
        for (uint64_t op : old_positions_chain)
        {
            old_hash_table->hash_table_array[op].unset_is_flagged_1();
        }
        // Then add k-1 previous k-mers so we can calculate hash values for all k-mers
        if (!break_out_because_no_previous)
        {
            //std::cout << "ADDING EXTRA K-MERS OT COMPLETE THE LAST ONE\n";
            for (int j = 0; j < kmer_len-1; j++)
            {
                if (old_hash_table->hash_table_array[chain_position].is_occupied())
                {
                    old_positions_chain.push_back(chain_position);
                    chain_length += 1;
                    //std::cout << "ALSO ADDED CHAIN POSITION " << chain_position << "\n";
                }
                else
                {
                    //std::cout << "COULD NOT ADD CHAIN POSITION " << chain_position << "\n";
                    break;
                }
                // If the previous k-mer exists, we add it to the list
                if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                {
                    chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                    //std::cout << "this one has a previous k-mer at " << chain_position << "\n";
                }
                // Otherwise break the loop early
                else
                {
                    //std::cout << "this one did not have previous k-mer\n";
                    //std::cout << "NO PREVIOUS K-MER EXISTS IN CAHIN POSITION " << chain_position << "\n";
                    if (old_hash_table->hash_table_array[chain_position].is_in_main_array())
                    {
                        //std::cout << "BUT IS IN MAIN ARRAY??????\n";
                    }
                    break;
                }
            }
        }
        

        // Now, start pushing characters from the chain (in reverse order) into the k-mer factory and calculate their new hash values
        
        // First, check if the last k-mer is in the temp array
        //std::cout << "FINDING FULL CHAIN for chain where initial length is: " << chain_length << "\n";
        uint64_t at_chain_position = chain_length-1;
        //std::cout << "THERE ARE THIS MANY ITEMS IN THE CHAIN VECTOR " << old_positions_chain.size() << "\n";
        uint64_t last_position_in_chain = old_positions_chain[at_chain_position];
        //std::cout << "The last position in chain is: " << last_position_in_chain << "\n";
        // If the last k-mer is in temp, load it fully into k-mer factory
        if (!old_hash_table->hash_table_array[last_position_in_chain].is_in_main_array())
        {
            //std::cout << "was not in main\n";
            // Get the slot in temp
            uint64_t old_temp_slot = old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot();
            for (int k = 0; k < kmer_factory->number_of_blocks; k++)
            {
                //std::cout << "POSSIBLE ERROR SPOT STARTS\n";
                uint64_t block_to_be_pushed = old_hash_table->temp_array[old_temp_slot*kmer_blocks+k];
                //std::cout << "GOT OVER IT\n";
                // If we are pushing the first block do this (possibly not full block)
                if (k == 0)
                {
                    //std::cout << "...\n";
                    for (int m = 0; m < kmer_factory->bits_in_last_block; m+=2)
                    {
                        //std::cout << "pushed char: " << (uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2))) << "\n"; 
                        kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2)));
                        new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                        block_to_be_pushed = block_to_be_pushed << 2;
                        
                    }
                    //std::cout << "...\n";
                }
                // Otherwise we are pushing a full block
                else 
                {
                    for (int m = 0; m < 64; m+=2)
                    {
                        kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (64-2)));
                        new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                        block_to_be_pushed = block_to_be_pushed << 2;
                    }
                }
            }
        } // Now the kmer_factory holds the "first" k-mer completely
        // Alternatively, if the last k-mer is in the main array do this
        else
        {
            //std::cout << "was in main\n";
            at_chain_position = chain_length-1;
            //std::cout << "...\n";
            for (int k = 0; k < kmer_len; k++)
            {
                //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                at_chain_position-=1;
            }
            //std::cout << "...\n";
            at_chain_position+=1;
        }

        //std::cout << "FULL CHAIN FOUND\n";
        //std::cout << "START ADDING FROM OLD POSITION " << old_positions_chain[at_chain_position] << "\n";

        // Now k-mer is loaded with a k-mer whose position in the old hash table is in chain_positions[at_chain_position]

        // Start checking the k-mers and add them to the new hash table it flag2 is not set
        uint64_t at_chain_position_backup = at_chain_position;
        while(at_chain_position >= 0)
        {
            // Get the position in the old array
            uint64_t old_position = old_positions_chain[at_chain_position];
            //std::cout << "- MIGRATING k-MER FROM OLD POSITION " << old_position << "\n";
            
            // If not flagged2, migrate to new hash table
            if (!old_hash_table->hash_table_array[old_position].is_flagged_2())
            {
                //std::cout << "-- the previous k-mer for this position lies at slot " << old_hash_table->hash_table_array[old_position].get_previous_kmer_slot() << "\n";
                //std::cout << "--- actually transfering position " << old_position << "\n";
                // Initial new position
                uint64_t new_position = new_hasher->get_current_hash();
                uint64_t probing_backup = new_position;
                // Probe for empty new position
                while(hash_table_array[new_position].is_occupied())
                {
                    new_position += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
                    if (new_position >= size)
                        new_position = new_position % size;
                    if (new_position == probing_backup){
                        std::cout << "RESIZING FAILURE, VERY LETHAL\n";
                        exit(1);
                    }
                }

                // Insert into the found empty position
                hash_table_array[new_position].set_data(old_hash_table->hash_table_array[old_position].get_data());
                hash_table_array[new_position].set_count(old_hash_table->hash_table_array[old_position].get_count());
                hash_table_array[new_position].set_previous_kmer_slot(old_position);
                hash_table_array[new_position].unset_is_flagged_1();
                hash_table_array[new_position].unset_is_flagged_2();
                hash_table_array[new_position].unset_is_written_in_output();
                
                // Remember where the old position was migrated to
                where_did_i_go[old_position] = new_position;
                // Set flag
                old_hash_table->hash_table_array[old_position].set_is_flagged_2();
                migrated_count+=1;

                // Count new slot
                // Add
                // Flags
            }
            else
            {
                //std::cout << "-- the previous k-mer for this position lies at slot " << hash_table_array[old_hash_table->hash_table_array[old_position].get_previous_kmer_slot()].get_previous_kmer_slot() << "\n";
                //std::cout << "--- not transfering position " << old_position <<  " because already transferred\n";
            }

            at_chain_position-=1;
            // All is done
            if (at_chain_position == -1)
            {
                break;
            }
            // Or continue with the next
            else
            {
                //std::cout << "...\n";
                //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                //std::cout << "...\n";
            }
        }
        //old_hash_table->hash_table_array[check_position].set_is_flagged_2();
        check_position+=1;
    }

    // SET THE CORRECT PREVIOUS K-MERS IN THE NEW HASH TABLE
    uint64_t check_position_old;
    uint64_t check_position_old_previous;
    uint64_t check_position_new_previous;
    check_position = 0;
    while (check_position < size)
    {
        if (hash_table_array[check_position].is_occupied())
        {
            check_position_old = hash_table_array[check_position].get_previous_kmer_slot();
            // If old k-mer was in main
            if (old_hash_table->hash_table_array[check_position_old].is_in_main_array())
            {
                check_position_old_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
                check_position_new_previous = where_did_i_go[check_position_old_previous];
            }
            // And if it was not
            else
            {
                check_position_new_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
            }
            
            hash_table_array[check_position].set_previous_kmer_slot(check_position_new_previous);
        }
        check_position += 1;
    }
    std::cout << "Resizing finished\n";
}

// ==============================================================================================================




// === For CANONICAL pointer hash table =========================================================================================


PointerHashTableCanonical::PointerHashTableCanonical(uint64_t s, uint64_t k, uint64_t b)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMer[size];
    bits_per_char = 2;
    inserted_items = 0;
    inserted_complete_kmers = 0;
    kmer_blocks = b;
    probe_hasher = new ProbeHasher1();
    probing_prime = mathfunctions::next_prime(uint64_t(std::floor(size/13.0)));
    max_temp_slots = 100;
    touched_temp_slots = 0;
    temp_slots_in_use = 0;
    max_temp_slot_in_use = 0;
    smallest_unused_temp_slot = 0;
    temp_array = std::vector<uint64_t>(b*max_temp_slots, uint64_t(0));
    temp_free_slots = std::vector<uint8_t>(max_temp_slots, 1);
}

PointerHashTableCanonical::~PointerHashTableCanonical()
{
    delete[] hash_table_array;
    delete probe_hasher;
}

uint64_t PointerHashTableCanonical::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(); 
}

bool PointerHashTableCanonical::kmer_in_slot_is_complete(uint64_t slot)
{
    return hash_table_array[slot].is_complete(); 
}

void PointerHashTableCanonical::resize()
{
    std::cout << "Hash table resizing not implemented yet...\n";
    exit(1);
}

uint64_t PointerHashTableCanonical::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t PointerHashTableCanonical::get_number_of_inserted_complete_kmers()
{
    return inserted_complete_kmers;
}

// Check if the currently loaded k-mer is in the given slot
// Return -1 if not in slot
// Return 0 if in slot but it is in the temp array from the very beginning
// Return 1 if in slot
int PointerHashTableCanonical::full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t slot)
{
    // Check if the slot is occupied
    if (!hash_table_array[slot].is_occupied())
        return -1;
    // Check if fully in temp
    if (!hash_table_array[slot].is_in_main_array())
    {
        for (int i = 0; i < kmer_factory->number_of_blocks; i++)
        {
            if (kmer_factory->get_canonical_block(i) != temp_array[hash_table_array[slot].get_previous_kmer_slot()*kmer_factory->number_of_blocks + i])
            {
                return -1;
            }  
        }
        return 0;
    }
    bool check_started_from_left = hash_table_array[slot]->check_starts_from_left(); 
    bool checking_left_side = check_started_from_left;
    bool flip_char = false;
    uint64_t previous_kmer_slot = 0;
    // If we pass the initial checks, we have to start comparing characters one by one
    uint64_t dbi = kmer_factory->number_of_blocks-1; // which data block we are looking at
    uint64_t db; // store current data block content here
    if (check_started_from_left)
    {
        db = kmer_factory->get_noncanonical_block(dbi); 
    }
    else
    {
        db = kmer_factory->get_canonical_block(dbi);
    }
    
    uint64_t cbc = 0; // number of checked data block characters
    uint64_t mc = 0; // matched characters
    bool check_rest_from_temp = false;
    uint64_t temp_array_slot = 0;
    uint64_t match_char;

    

    while (mc < kmer_factory->get_number_of_stored_characters())
    {
        // If the current hash table slot is unoccupied the k-mer does not exist
        if (!hash_table_array[slot].is_occupied())
        {
            return -1;
        }
        if (!hash_table_array[slot].is_in_main_array())
        {
            //std::cout << "Found link to temp mid search, start hybrid search\n";
            check_rest_from_temp = true;
            temp_array_slot = hash_table_array[slot].get_previous_kmer_slot();
            //std::cout << "Hybrid check starts in slot " << temp_array_slot << "\n";
            break;
        }
        // Check that the characters match
        match_char = (db&uint64_t(3));
        if (check_started_from_left)
        {
            match_char = uint64_t(3) - match_char;
        }
        if (uint64_t(hash_table_array[slot].get_character(checking_left_side, flip_char)) != match_char)
        {
            //std::cout << "k-mer not in this slot, normal check\n";
            return -1;
        // Update to look at next characters if the current ones match
        } 
        else 
        {
            mc+=1;
            if (mc == kmer_factory->get_number_of_stored_characters())
                return 1;
            if (hash_table_array[slot].prev_kmer_exists()){
                previous_kmer_slot = hash_table_array[slot].get_previous_kmer_slot();
                // Change character side if different orientation
                if (!hash_table_array[slot]->same_orientation_as_previous())
                {
                    checking_left_side = !checking_left_side;
                }
                // Change character flipping if necessary
                flip_char = false;
                if (check_started_from_left == hash_table_array[slot]->check_starts_from_left())
                {
                    if (!hash_table_array[slot]->same_orientation_as_previous())
                    {
                        flip_char = true;
                    }
                } 
                else 
                {
                    if (hash_table_array[slot]->same_orientation_as_previous())
                    {
                        flip_char = true;
                    }
                }
                // Character side and flipping resolved, keep going normally
                slot = previous_kmer_slot;
            } else {
                //std::cout << "VERY SUS ERROR!! Should not be possible to get here\n";
                return -1;
            }
            cbc += 1;
            db >>= 2;
            if (cbc >= 32){
                cbc = 0;
                dbi -= 1;
                //db = kmer_factory->get_canonical_block(dbi);
                if (check_started_from_left)
                {
                    db = kmer_factory->get_noncanonical_block(dbi); 
                }
                else
                {
                    db = kmer_factory->get_canonical_block(dbi);
                }
            }
        }
    }
    // If the rest of he characters must be checked from temp
    if (check_rest_from_temp)
    {
        if (checking_left_side)
        {
            uint64_t current_temp_array_position = temp_array_slot*kmer_factory->number_of_blocks;
            uint64_t temp_array_block = temp_array[current_temp_array_position];
            int unchecked_temp_block_bits = kmer_factory->bits_in_last_block;
            temp_array_block <<= (64-unchecked_temp_block_bits);
            //std::cout << "Hybrid check adjusted slot is " << current_temp_array_position << "\n";
            //std::cout << "and the block in question is " << temp_array_block << "\n";
            //std::cout << "and the kmer factory block is now " << db << "\n";
            while (mc < kmer_factory->get_number_of_stored_characters())
            {
                //std::cout << "COMPARING FACTORY " << (db&uint64_t(3)) << " WITH TEMP " << (temp_array_block&uint64_t(3)) << "\n";
                match_char = (db&uint64_t(3));
                if (check_started_from_left)
                {
                    match_char = uint64_t(3) - match_char;
                }
                if (flip_char)
                {
                    if (match_char != (uint64_t(3)-((temp_array_block>>62)&uint64_t(3))))
                    {
                        //std::cout << "k-mer not in this slot, hybrid check\n";
                        return -1;
                    }
                }
                else
                {
                    if (match_char != ((temp_array_block>>62)&uint64_t(3)))
                    {
                        //std::cout << "k-mer not in this slot, hybrid check\n";
                        return -1;
                    }
                }
                                
                mc+=1;
                if (mc == kmer_factory->get_number_of_stored_characters())
                {
                    //std::cout << "Found using PARTIAL TEMP\n";
                    //std::cout << "Search hit during hybrid search\n";
                    return 1;
                }
                    
                cbc+=1;
                db>>=2;
                if (cbc >= 32){
                    cbc = 0;
                    dbi -= 1;
                    //db = kmer_factory->get_canonical_block(dbi);
                    if (check_started_from_left)
                    {
                        db = kmer_factory->get_noncanonical_block(dbi); 
                    }
                    else
                    {
                        db = kmer_factory->get_canonical_block(dbi);
                    }
                }
                temp_array_block <<= 2;
                unchecked_temp_block_bits-=2;
                if (unchecked_temp_block_bits == 0)
                {
                    unchecked_temp_block_bits = 64;
                    current_temp_array_position += 1;
                    temp_array_block = temp_array[current_temp_array_position];
                }
            }
        }
        else
        {
            uint64_t current_temp_array_position = temp_array_slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1;
            uint64_t temp_array_block = temp_array[current_temp_array_position];
            int unchecked_temp_block_bits = 64;
            //std::cout << "Hybrid check adjusted slot is " << current_temp_array_position << "\n";
            //std::cout << "and the block in question is " << temp_array_block << "\n";
            //std::cout << "and the kmer factory block is now " << db << "\n";
            while (mc < kmer_factory->get_number_of_stored_characters())
            {
                //std::cout << "COMPARING FACTORY " << (db&uint64_t(3)) << " WITH TEMP " << (temp_array_block&uint64_t(3)) << "\n";
                match_char = (db&uint64_t(3));
                if (check_started_from_left)
                {
                    match_char = uint64_t(3) - match_char;
                }
                if (flip_char)
                {
                    if (match_char != (uint64_t(3)-(temp_array_block&uint64_t(3))))
                    {
                        //std::cout << "k-mer not in this slot, hybrid check\n";
                        return -1;
                    }
                }
                else
                {
                    if (match_char != (temp_array_block&uint64_t(3)))
                    {
                        //std::cout << "k-mer not in this slot, hybrid check\n";
                        return -1;
                    }
                }
                mc+=1;
                if (mc == kmer_factory->get_number_of_stored_characters())
                {
                    //std::cout << "Found using PARTIAL TEMP\n";
                    //std::cout << "Search hit during hybrid search\n";
                    return 1;
                }
                    
                cbc+=1;
                db>>=2;
                if (cbc >= 32){
                    cbc = 0;
                    dbi -= 1;
                    //db = kmer_factory->get_canonical_block(dbi);
                    if (check_started_from_left)
                    {
                        db = kmer_factory->get_noncanonical_block(dbi); 
                    }
                    else
                    {
                        db = kmer_factory->get_canonical_block(dbi);
                    }
                }
                temp_array_block >>= 2;
                unchecked_temp_block_bits-=2;
                if (unchecked_temp_block_bits == 0)
                {
                    unchecked_temp_block_bits = 64;
                    current_temp_array_position -= 1;
                    temp_array_block = temp_array[current_temp_array_position];
                }
            }
        }
    }
    // If the characters have matched until now, the k-mer exists in the hash table
    return 1;
}



// FULL K-MER SLOT CHECK MODIFIED, CONTINUE MODIFICATIONS FROM HERE




// Find where a k-mer is in the hash table
std::tuple<bool,uint64_t> PointerHashTableCanonical::find_and_modify_using_prev(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = kmer_hash;
    // Probe as long as an occupied slot is found
    while (hash_table_array[kmer_slot].is_occupied())
    {
        // If previous k-mer exists, make a quick previous k-mer check
        if (previous_kmer_exists)
        {
            if (uint64_t(hash_table_array[kmer_slot].get_canonical_right_character() == uint64_t(kmer_factory->get_canonical_right_character())))
            {
                if (previous_kmer_slot == hash_table_array[kmer_slot].get_previous_kmer_slot())
                {
                    return std::make_tuple(true,kmer_slot);
                }
            }
        }
        // If quick check cannot be done or it is not successful, do a full check
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        // If k-mer in slot and the slot does not point to temp immediately
        if (check_result == 1)
        {
            return std::make_tuple(false,kmer_slot);
        }
        // If k-mer is not in the curent slot, probe to the next slot
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe_1(kmer_factory->get_canonical_block(kmer_factory->number_of_blocks-1), probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                //resize();
                exit(1);
            }
        }
        // If k-mer is in the slot, but the slot points to temp
        else if (check_result == 0)
        {
            // If we have a previous k-mer, migrate from temp to main
            if (previous_kmer_exists)
            {
                uint64_t temp_array_slot_in_question = hash_table_array[kmer_slot].get_previous_kmer_slot();
                // Release temp slot and recalculate
                temp_free_slots[temp_array_slot_in_question] = 1;
                temp_slots_in_use -= 1;
                for (int i = 0; i < kmer_factory->number_of_blocks; i++)
                {
                    temp_array[temp_array_slot_in_question*kmer_factory->number_of_blocks + i] = 0;
                }
                // Add main array count
                inserted_items+=1;
                hash_table_array[kmer_slot].set_prev_kmer_exists();
                hash_table_array[kmer_slot].set_previous_kmer_slot(previous_kmer_slot);   
                // Increase count
                //hash_table_array[kmer_slot].increase_count();
                // Set last character
                //hash_table_array[kmer_slot].set_character(kmer_factory->get_newest_character());
                hash_table_array[kmer_slot].set_right_character(kmer_factory->get_canonical_right_character());
                hash_table_array[kmer_slot].set_left_character(kmer_factory->get_canonical_left_character());
                // Set occupied
                hash_table_array[kmer_slot].set_occupied();
                // Set in main array
                hash_table_array[kmer_slot].set_in_main_array();
                // Set complete/partial
                hash_table_array[kmer_slot].set_complete();
                inserted_complete_kmers += 1;
                // Set orientation
                if !(kmer_factory->previous_kmer_exists)
                {
                    std::cout << "HORRIBLE ERROR, SHOULD NOT BE POSSIBLE TO GET HERE\n";
                    exit(1);
                }
                if (kmer_factory->forward_is_canonical)
                {
                    hash_table_array[kmer_slot].unset_check_starts_from_left();
                    if (kmer_factory->previous_forward_was_canonical)
                    {
                        hash_table_array[kmer_slot].set_same_orientation_as_previous();
                    }
                    else
                    {
                        hash_table_array[kmer_slot].unset_same_orientation_as_previous();
                    }
                }
                else
                {
                    hash_table_array[kmer_slot].set_check_starts_from_left();
                    if (kmer_factory->previous_forward_was_canonical)
                    {
                        hash_table_array[kmer_slot].unset_same_orientation_as_previous();
                    }
                    else
                    {
                        hash_table_array[kmer_slot].set_same_orientation_as_previous();
                    }
                }
                return std::make_tuple(false,kmer_slot);
            }
            // If the new k-mer does not have a previous k-mer, do not migrate and instead increase count in temp
            else
            {
                return std::make_tuple(false,kmer_slot);
            }
        }
        else
        {
            std::cout << "ERROR : It should not be possible to end up in this situation\n";
            exit(1);
        }    
    }
    //std::cout << "SLOT " << kmer_slot <<  " WAS UNOCCUPIED\n";
    return std::make_tuple(false,size);
}

// Find where a k-mer is in the hash table
uint64_t PointerHashTableCanonical::find(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    while (hash_table_array[kmer_slot].is_occupied())
    {
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        if (check_result == 1)
        {
            return kmer_slot;
        } 
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe_1(kmer_factory->get_canonical_block(kmer_factory->number_of_blocks-1), probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                resize();
                exit(1);
            }
        }    
        else if (check_result == 0)
        {
            return kmer_slot;
        }
        else
        {
            std::cout << "ERROR : It should not be possible to end up in this situation\n";
            exit(1);
        }    
    }
    //std::cout << "SLOT " << kmer_slot <<  " WAS UNOCCUPIED\n";
    return size;
}

// Insert a new k-me to the hash table, we assume the k-mer is not already in the hash table
uint64_t PointerHashTableCanonical::insert_new_kmer_in_main(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, uint64_t prev_kmer_slot)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    
    while(hash_table_array[kmer_slot].is_occupied())
    {
        kmer_slot += probe_hasher->probe_1(kmer_factory->get_canonical_block(kmer_factory->number_of_blocks-1), probing_prime);
        //kmer_slot += 1;
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == kmer_hash)
        {
            //std::cout << "Hash table is full, needs to be resized\n";
            resize();
            exit(1);
        }
    }

    // Next, add the k-mer info to the hash table

    // First, determine some characteristics of the k-mer
    bool complete = kmer_factory->get_number_of_stored_characters() >= kmer_len;
    bool prev_exists = kmer_factory->get_number_of_stored_characters() > 1;

    inserted_items+=1;

    // Set previous k-mer status (and previous k-mer slot) 
    if (prev_exists){
        hash_table_array[kmer_slot].set_prev_kmer_exists();
        hash_table_array[kmer_slot].set_previous_kmer_slot(prev_kmer_slot);
    } else {
        //std::cout << "DISASTROUS ERROR\n\n\n\n";
        hash_table_array[kmer_slot].unset_prev_kmer_exists();
    }
    // Set count to one
    hash_table_array[kmer_slot].set_count(1);
    // Set last character
    //hash_table_array[kmer_slot].set_character(kmer_factory->get_newest_character());
    hash_table_array[kmer_slot].set_right_character(kmer_factory->get_canonical_right_character());
    hash_table_array[kmer_slot].set_left_character(kmer_factory->get_canonical_left_character());
    if (kmer_factory->forward_is_canonical)
    {
        hash_table_array[kmer_slot].unset_check_starts_from_left();
        if (kmer_factory->previous_forward_was_canonical)
        {
            hash_table_array[kmer_slot].set_same_orientation_as_previous();
        }
        else
        {
            hash_table_array[kmer_slot].unset_same_orientation_as_previous();
        }
    }
    else
    {
        hash_table_array[kmer_slot].set_check_starts_from_left();
        if (kmer_factory->previous_forward_was_canonical)
        {
            hash_table_array[kmer_slot].unset_same_orientation_as_previous();
        }
        else
        {
            hash_table_array[kmer_slot].set_same_orientation_as_previous();
        }
    }
    // Set occupied
    hash_table_array[kmer_slot].set_occupied();
    // Set in main array
    hash_table_array[kmer_slot].set_in_main_array();
    // Set complete/partial
    if (complete){
        hash_table_array[kmer_slot].set_complete();
        inserted_complete_kmers += 1;
    } else {
        hash_table_array[kmer_slot].unset_complete();
    }
    //std::cout << "Inserted " << hash_table_array[kmer_slot].get_data() << "\n";
    //std::cout << "INSERTED IN MAIN USING SLOT " << kmer_slot << " IN MAIN AND PREV LINK TO SLOT " << prev_kmer_slot << "\n";
    return kmer_slot;
}

// Insert a new k-me to the hash table, we assume the k-mer is not already in the hash table
uint64_t PointerHashTableCanonical::insert_new_kmer_in_temp(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    
    while(hash_table_array[kmer_slot].is_occupied())
    {
        kmer_slot += probe_hasher->probe_1(kmer_factory->get_canonical_block(kmer_factory->number_of_blocks-1), probing_prime);
        //kmer_slot += 1;
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == kmer_hash)
        {
            std::cout << "Hash table is full, needs to be resized\n";
            resize();
            exit(1);
        }
    }

    // Find the next free slot
    bool empty_temp_slot_found = false;
    for (int j = 0; j < max_temp_slots; j++)
    {
        if (temp_free_slots[j] == 1)
        {
            //std::cout << "Could be improved...\n";
            smallest_unused_temp_slot = j;
            empty_temp_slot_found = true;
            break;
        }
    }
    if (!empty_temp_slot_found)
    {
        std::cout << "TEMO SLOT RESIZING ERROR THAT SHOULD NOT HAPPEN\n";
        exit(1);
    }
    max_temp_slot_in_use = std::max(max_temp_slot_in_use, smallest_unused_temp_slot+1);
    temp_free_slots[smallest_unused_temp_slot] = 0;
    temp_slots_in_use += 1;
    
    hash_table_array[kmer_slot].set_occupied();
    hash_table_array[kmer_slot].unset_prev_kmer_exists();
    hash_table_array[kmer_slot].set_count(1);
    hash_table_array[kmer_slot].set_complete();
    hash_table_array[kmer_slot].unset_in_main_array();
    hash_table_array[kmer_slot].set_previous_kmer_slot(smallest_unused_temp_slot);
    for (int i = 0; i < kmer_factory->number_of_blocks; i++)
    {
        temp_array[smallest_unused_temp_slot*kmer_factory->number_of_blocks + i] = kmer_factory->get_canonical_block(i);
        //std::cout << "WE PUT IN TEMP THE FOLLOWING K-MER: " << kmer_factory->blocks[i] << "\n";
    }
    if (kmer_factory->forward_is_canonical)
    {
        hash_table_array[kmer_slot].unset_check_starts_from_left();
    }
    else
    {
        hash_table_array[kmer_slot].set_check_starts_from_left();
    }
    if (kmer_factory->previous_kmer_exists)
    {
        std::cout  << "ERROR 345, should not be able to get in here\n";
        exit(1);
    }

    if (smallest_unused_temp_slot == max_temp_slots-1)
    {
        std::cout << "Temp array resizing needed\n";
        std::cout << "Might be BUGGED.............................................\n";
        //exit(1);
        max_temp_slots *= 2;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
        //temp_free_slots = sdsl::bit_vector(max_temp_slots, 1);
        //temp_occupied_vector.resize(max_temp_slots);

    }
    //std::cout << "INSERTED IN TEMP USING SLOT " << kmer_slot << " IN MAIN AND SLOT " << used_temp_slots-1 << " IN TEMP\n";
    return kmer_slot;
}

uint64_t PointerHashTableCanonical::find_and_increment(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_ip_step)
{
    bool using_prev;
    uint64_t kmer_slot;
    if (clean_ip_step)
    {
        kmer_slot = find(kmer_factory, kmer_hash);
    }
    else
    {
        std::tie(using_prev, kmer_slot) = find_and_modify_using_prev(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot);
    }
    if (kmer_slot < size)
    {
        if (kmer_factory->get_number_of_stored_characters() == kmer_len)
        {
            hash_table_array[kmer_slot].increase_count();
        }
        else
        {
            //std::cout << "===== SUS INCREASE TRY IN 'kmer_hash_table.cpp' function 'find_and_increment' =====\n";
        }
        
        if(previous_kmer_exists)
        {   
            if (hash_table_array[hash_table_array[kmer_slot].get_previous_kmer_slot()].get_count() < hash_table_array[previous_kmer_slot].get_count())
            {
                //std::cout << "PREVIOUS K-MER SWAPPED\n";
                hash_table_array[kmer_slot].set_previous_kmer_slot(previous_kmer_slot);
            }
        }
    }
        
    return kmer_slot;
}


uint64_t PointerHashTableCanonical::add_kmer(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot, bool clean_up_step)
{
    // Try to find and increment
    uint64_t kmer_slot = find_and_increment(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot, clean_up_step);
    // If that did not work, insert new
    if (kmer_slot == size)
    {
        // If previous k-mer exists, insert normally
        if (previous_kmer_exists || clean_up_step)
        {
            kmer_slot = insert_new_kmer_in_main(kmer_factory, kmer_hash, previous_kmer_slot);
            //std::cout << "Inserting new to main array\n";
        }
        // If previous k-mer does not exist, insert in temp
        else
        {
            kmer_slot = insert_new_kmer_in_temp(kmer_factory, kmer_hash);
            //std::cout << "Inserting new to temp array\n";
        }
    }
    else
    {
        //std::cout << "k-mer was found and its count increased\n";
    }
    return kmer_slot;
}



/*






Continue rechecking modifications here
    - Modification needed because: check starting from left must be checked in reverse character order












*/







void PointerHashTableCanonical::clean_up_temp_array(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher)
{
    uint64_t brought_back_home = 0;
    for (int i = 0; i < size; i++)
    {
        if (!hash_table_array[i].is_occupied())
        {
            continue;
        }
        if (hash_table_array[i].is_in_main_array())
        {
            continue;
        }
        brought_back_home+=1;
        //std::cout << "K-MERS BROUGHT BACK HOME " << brought_back_home << "\n";
        //std::cout << "SLOT " << i << " IN MAIN WAS STORED IN TEMP, BRINGING IT BACK HOME\n";
        temp_slots_in_use -= 1;
        kmer_factory->reset();
        hasher->reset();

        uint64_t handled_characters = 0;
        uint64_t previous_kmer_slot = size;
        bool previous_kmer_exists = false;
        uint64_t slot_for_me = 0;
        uint64_t current_slot_in_temp = hash_table_array[i].get_previous_kmer_slot();

        while (handled_characters < kmer_len-1)
        {
            if (handled_characters > 0)
            {
                previous_kmer_exists = true;
            }
            kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
            hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
            if (kmer_factory->forward_kmer_is_canonical())
            {
                slot_for_me = hasher->get_current_hash_forward();
            }
            else
            {
                slot_for_me = hasher->get_current_hash_backward();
            }
            slot_for_me = find_and_increment_andifnot_insert(kmer_factory, slot_for_me, previous_kmer_exists, previous_kmer_slot, true);
            previous_kmer_slot = slot_for_me;
            handled_characters += 1;
        }
        slot_for_me = i;
        kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
        hasher->update_rolling_hash(kmer_factory->get_forward_newest_character(), kmer_factory->get_forward_pushed_off_character());
        if (hasher->get_current_hash() != i)
        {
            //std::cout << "SUS ERROR !!\n";
        }

        hash_table_array[slot_for_me].set_prev_kmer_exists();
        hash_table_array[slot_for_me].set_previous_kmer_slot(previous_kmer_slot);
        //hash_table_array[slot_for_me].set_character(kmer_factory->get_newest_character());
        hash_table_array[kmer_slot].set_right_character(kmer_factory->get_canonical_right_character());
        hash_table_array[kmer_slot].set_left_character(kmer_factory->get_canonical_left_character());
        if (kmer_factory->forward_is_canonical)
        {
            hash_table_array[kmer_slot].unset_check_starts_from_left();
            if (kmer_factory->previous_forward_was_canonical)
            {
                hash_table_array[kmer_slot].set_same_orientation_as_previous();
            }
            else
            {
                hash_table_array[kmer_slot].unset_same_orientation_as_previous();
            }
        }
        else
        {
            hash_table_array[kmer_slot].set_check_starts_from_left();
            if (kmer_factory->previous_forward_was_canonical)
            {
                hash_table_array[kmer_slot].unset_same_orientation_as_previous();
            }
            else
            {
                hash_table_array[kmer_slot].set_same_orientation_as_previous();
            }
        }
        hash_table_array[slot_for_me].set_in_main_array();
        hash_table_array[slot_for_me].set_complete();
        inserted_complete_kmers += 1;
        inserted_items+=1;
    }
}


char PointerHashTableCanonical::get_oldest_char_in_temp_slot_and_shift_left(uint64_t slot, KMerFactoryCanonical2BC* kmer_factory)
{
    uint64_t next_char_chunk = ((temp_array[slot*kmer_factory->number_of_blocks]) >> (kmer_factory->bits_in_last_block-2));
    char next_char = twobitstringfunctions::int2char(next_char_chunk);
    for (int i = 0; i < kmer_factory->number_of_blocks-1; i++)
    {
        temp_array[slot*kmer_factory->number_of_blocks + i] = temp_array[slot*kmer_factory->number_of_blocks + i] << 2;
        temp_array[slot*kmer_factory->number_of_blocks + i] = (temp_array[slot*kmer_factory->number_of_blocks + i] | (temp_array[slot*kmer_factory->number_of_blocks + i + 1] >> 64 - 2));
    }
    temp_array[slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1] = (temp_array[slot*kmer_factory->number_of_blocks + kmer_factory->number_of_blocks-1] << 2);
    temp_array[slot*kmer_factory->number_of_blocks] = (temp_array[slot*kmer_factory->number_of_blocks] & (kmer_factory->used_left_block_mask));
    return next_char;
}

uint64_t PointerHashTableCanonical::get_number_of_max_temp_slots()
{
    return max_temp_slots;
}

uint64_t PointerHashTableCanonical::get_number_of_touched_temp_slots()
{
    return touched_temp_slots;
}

uint64_t PointerHashTableCanonical::get_number_of_temp_slots_in_use()
{
    return temp_slots_in_use;
}

uint64_t PointerHashTableCanonical::get_number_of_max_temp_slots_in_use()
{
    return max_temp_slot_in_use;
}

uint64_t PointerHashTableCanonical::get_smallest_unused_temp_slot()
{
    for (int j = 0; j < max_temp_slots; j++)
    {
        if (temp_free_slots[j] == 1)
        {
            smallest_unused_temp_slot = j;
            break;
        }
    }
    return smallest_unused_temp_slot;
}

void PointerHashTableCanonical::write_kmers_on_disk(KMerFactory2BC* kmer_factory, uint64_t min_abundance, bool only_canonical, std::string output_path)
{
    //std::string output_file_path = output_directory + "/" + std::to_string(kmer_len) + "-mers_" + std::to_string(min_abundance) + "-abundance_2-mode.txt";
	std::ofstream output_file(output_path);

    uint64_t check_position = 0;
    uint64_t chain_position = 0;
    std::vector<uint64_t> current_chain_abundances;
    //std::vector<bool> current_chain_completeness;
    std::string current_chain_string = "";
    int current_chain_streak = 0;
    int round_one_kmers = 0;
    int round_two_kmers = 0;

    sdsl::bit_vector is_referenced(size,0);
    while(check_position < size)
    {
        if (hash_table_array[check_position].is_occupied())
        {
            is_referenced[hash_table_array[check_position].get_previous_kmer_slot()] = 1;
        }
        check_position+=1;
    }

    for (int round = 0; round < 2; round++)
    {
        //if (round == 0)
        //{
        //    continue;
        //}          
        check_position = 0;
        chain_position = 0;
        current_chain_abundances.clear();
        current_chain_string = "";
        current_chain_streak = 0;

        while (check_position < size)
        {
            if (round == 0)
            {
                if (is_referenced[check_position] == 1)
                {
                    check_position+=1;
                    continue;
                }
            }
            
            // Written in output == is flagged_1
            if ((hash_table_array[check_position].is_complete()) && (!hash_table_array[check_position].is_flagged_1()) && (hash_table_array[check_position].is_occupied()) && (hash_table_array[check_position].is_complete()) && (hash_table_array[check_position].get_count() >= min_abundance))
            {
                bool chain_check_starts_from_left = hash_table_array[check_position].check_starts_from_left();
                bool current_char_taken_from_left = chain_check_starts_from_left;
                bool flip_char = false;
                bool can_continue_chain = true;

                chain_position = check_position;
                // Chain together all k-mers that have not been written into output
                while (!hash_table_array[chain_position].is_flagged_1() && (hash_table_array[chain_position].is_complete()))
                {
                    current_chain_streak+=1;
                    current_chain_abundances.push_back(hash_table_array[chain_position].get_count());
                    //current_chain_completeness.push_back(hash_table_array[chain_position].is_complete());
                    current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character(current_char_taken_from_left, flip_char))) + current_chain_string;
                    //current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character())) + current_chain_string;
                    hash_table_array[chain_position].set_is_flagged_1();
                    if(hash_table_array[chain_position].prev_kmer_exists())
                    {
                        // Change character side if different orientation
                        if (!hash_table_array[chain_position]->same_orientation_as_previous())
                        {
                            current_char_taken_from_left = !current_char_taken_from_left;
                        }
                        // Change character flipping if necessary
                        flip_char = false;
                        if (chain_check_starts_from_left == hash_table_array[chain_position]->check_starts_from_left())
                        {
                            if (!hash_table_array[chain_position]->same_orientation_as_previous())
                            {
                                flip_char = true;
                            }
                        } 
                        else 
                        {
                            if (hash_table_array[chain_position]->same_orientation_as_previous())
                            {
                                flip_char = true;
                            }
                        }
                        chain_position = hash_table_array[chain_position].get_previous_kmer_slot();
                    }
                    else
                    {
                        can_continue_chain = false;
                        break;
                    }
                }
                // Complete the last k-mer in the chain
                if (can_continue_chain)
                {
                    for (int i = 0; i < kmer_len-1; i++)
                    {
                        //current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character())) + current_chain_string;
                        current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character(current_char_taken_from_left, flip_char))) + current_chain_string;
                        if (hash_table_array[chain_position].prev_kmer_exists())
                        {
                            // Change character side if different orientation
                            if (!hash_table_array[chain_position]->same_orientation_as_previous())
                            {
                                current_char_taken_from_left = !current_char_taken_from_left;
                            }
                            // Change character flipping if necessary
                            flip_char = false;
                            if (chain_check_starts_from_left == hash_table_array[chain_position]->check_starts_from_left())
                            {
                                if (!hash_table_array[chain_position]->same_orientation_as_previous())
                                {
                                    flip_char = true;
                                }
                            } 
                            else 
                            {
                                if (hash_table_array[chain_position]->same_orientation_as_previous())
                                {
                                    flip_char = true;
                                }
                            }
                            chain_position = hash_table_array[chain_position].get_previous_kmer_slot();
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                
                // Write the chained k-mer into output
                for (int j = 0; j < current_chain_streak; j++)
                {
                    if (j+kmer_len > current_chain_string.length())
                    {
                        break;
                    }
                    if (round == 0)
                        round_one_kmers+=1;
                    else
                        round_two_kmers+=1;

                    if ((current_chain_abundances[current_chain_streak-j-1] >= min_abundance))
                    {
                        std::string current_kmer = current_chain_string.substr(j, kmer_len);
                        if (!purestringfunctions::is_canonical(current_kmer))
                        {
                            purestringfunctions::reverse_this_string(current_kmer);
                        }
                        output_file << current_kmer << " " << std::to_string(current_chain_abundances[current_chain_streak-j-1]) << "\n";
                    }
                }
                // Reset 
                //current_chain_completeness.clear();
                current_chain_abundances.clear();
                current_chain_string = "";
                current_chain_streak = 0;
            }
            check_position += 1;
        }
    }
    std::cout << "k-mer starts in round one: " << round_one_kmers << "\n";
    std::cout << "k-mer starts in round two: " << round_two_kmers << "\n";

    output_file.close(); // close file
	output_file.clear(); // clear flags
} 

uint64_t PointerHashTableCanonical::get_size()
{
    return size;
}


// Last thing to modify !!!!!!!!
void PointerHashTableCanonical::copy_content_from_another_hash_table(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* new_hasher, PointerHashTableCanonical* old_hash_table)
{
    if (max_temp_slots != old_hash_table->max_temp_slots)
    {
        max_temp_slots = old_hash_table->max_temp_slots;
        temp_array.resize(kmer_blocks*max_temp_slots, 0);
        temp_free_slots.resize(max_temp_slots, 1);
    }
    // Copy content from temp array
    for (int i = 0; i < old_hash_table->temp_array.size(); i++)
    {
        this->temp_array[i] = old_hash_table->temp_array[i];
    }
    // Also remember to copy the bit vector here
    for (int i = 0; i < old_hash_table->temp_free_slots.size(); i++)
    {
        this->temp_free_slots[i] = old_hash_table->temp_free_slots[i];
    }
    inserted_items = old_hash_table->inserted_items;
    inserted_complete_kmers = old_hash_table->inserted_complete_kmers;
    touched_temp_slots = old_hash_table->touched_temp_slots;
    temp_slots_in_use = old_hash_table->temp_slots_in_use;
    smallest_unused_temp_slot = old_hash_table->smallest_unused_temp_slot;

    

    // Start checking slots in the old hash table
    uint64_t check_position = 0;
    uint64_t chain_position;
    uint64_t chain_length;
    std::vector<uint64_t> old_positions_chain;
    //std::vector<uint64_t> where_did_i_go(old_hash_table->size, 0);
    int percentage_checked = 0;
    uint64_t migrated_count = 0;
    
    sdsl::bit_vector is_referenced(old_hash_table->get_size(),0);
    while(check_position < old_hash_table->get_size())
    {
        if (old_hash_table->hash_table_array[check_position].is_occupied())
        {
            is_referenced[old_hash_table->hash_table_array[check_position].get_previous_kmer_slot()] = 1;
        }
        check_position+=1;
    }

    int round_one_kmers = 0;
    int round_two_kmers = 0;

    // Flag 1 used in loop 
    // FLag 2 used to check if already transferred
    int percentage = 1;

    for (int round = 0; round < 2; round++)
    {
        check_position = 0;

        while (check_position < old_hash_table->get_size())
        {
            if (round == 0)
            {
                if (is_referenced[check_position] == 1)
                {
                    check_position+=1;
                    continue;
                }
            }

            //std::cout << "CHECK POSITION " << check_position << "\n";
            if (((100.0*migrated_count) / old_hash_table->get_size()) >=percentage)
            {
                std::cout << "Resizing... " << percentage << "%" << " done\n";
                percentage+=1;
            }
            //if (check_position % 100000 == 0)
            //    std::cout << "Resizing... " << check_position << " / " << old_hash_table->size << "\n";
            // If the k-mer at the check position is already migrated, continue to next position
            if (old_hash_table->hash_table_array[check_position].is_flagged_2())
            {
                //std::cout << "## k-mer already migrated at position " << check_position << "\n";
                check_position += 1;
                continue;
            }
            if (!old_hash_table->hash_table_array[check_position].is_occupied())
            {
                //std::cout << "?? no k-mer at position " << check_position << "\n";
                old_hash_table->hash_table_array[check_position].set_is_flagged_2();
                check_position += 1;
                migrated_count+=1;
                continue;
            }
            // Reset the previous chain
            old_positions_chain.clear();
            kmer_factory->reset();
            new_hasher->reset();

            // If the k-mer has not yet been migrated, create the chain and migrate all
            chain_position = check_position;
            chain_length = 0;
            // Chain all k-mers that have not been migrated AND that are not laready in the loop
            bool break_out_because_no_previous = false;
            while((!old_hash_table->hash_table_array[chain_position].is_flagged_1()) && (!old_hash_table->hash_table_array[chain_position].is_flagged_2()))
            {
                //std::cout << "CHAIN POSITION IS " << chain_position << "\n";
                //std::cout << "- its previous k-mer is at position " << old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot() << "\n";
                // Add the chain position to the chain position list
                old_positions_chain.push_back(chain_position);
                chain_length += 1;
                // Mark the k-mer so that we do not include it twice in the loop
                old_hash_table->hash_table_array[chain_position].set_is_flagged_1();
                // If previous k-mer exists, we check it next
                if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                {
                    if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                    {
                        chain_position = hash_table_array[old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot()].get_previous_kmer_slot();
                        std::cout << "SUS BEHAVIOR IF WE GET HERE\n";
                    }
                    else
                    {
                        chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                    }
                }
                // Otherwise we leave the loop
                else
                {
                    //std::cout << "~~ this position does not have previous k-mer\n";
                    break_out_because_no_previous = true;
                    break;
                }
            }
            // Unset loop flags
            for (uint64_t op : old_positions_chain)
            {
                old_hash_table->hash_table_array[op].unset_is_flagged_1();
            }
            // Then add k-1 previous k-mers so we can calculate hash values for all k-mers
            if (!break_out_because_no_previous)
            {
                //std::cout << "ADDING EXTRA K-MERS OT COMPLETE THE LAST ONE\n";
                for (int j = 0; j < kmer_len-1; j++)
                {
                    if (old_hash_table->hash_table_array[chain_position].is_occupied())
                    {
                        old_positions_chain.push_back(chain_position);
                        chain_length += 1;
                        //std::cout << "ALSO ADDED CHAIN POSITION " << chain_position << "\n";
                    }
                    else
                    {
                        //std::cout << "COULD NOT ADD CHAIN POSITION " << chain_position << "\n";
                        break;
                    }
                    // If the previous k-mer exists, we add it to the list
                    if (old_hash_table->hash_table_array[chain_position].prev_kmer_exists())
                    {
                        if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                        {
                            chain_position = hash_table_array[old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot()].get_previous_kmer_slot();
                        }
                        else
                        {
                            chain_position = old_hash_table->hash_table_array[chain_position].get_previous_kmer_slot();
                        }
                        //std::cout << "this one has a previous k-mer at " << chain_position << "\n";
                    }
                    // Otherwise break the loop early
                    else
                    {
                        //std::cout << "this one did not have previous k-mer\n";
                        //std::cout << "NO PREVIOUS K-MER EXISTS IN CAHIN POSITION " << chain_position << "\n";
                        if (old_hash_table->hash_table_array[chain_position].is_in_main_array())
                        {
                            //std::cout << "BUT IS IN MAIN ARRAY??????\n";
                        }
                        break;
                    }
                }
            }
            

            // Now, start pushing characters from the chain (in reverse order) into the k-mer factory and calculate their new hash values
            
            // First, check if the last k-mer is in the temp array
            //std::cout << "FINDING FULL CHAIN for chain where initial length is: " << chain_length << "\n";
            uint64_t at_chain_position = chain_length-1;
            //std::cout << "THERE ARE THIS MANY ITEMS IN THE CHAIN VECTOR " << old_positions_chain.size() << "\n";
            uint64_t last_position_in_chain = old_positions_chain[at_chain_position];
            //std::cout << "The last position in chain is: " << last_position_in_chain << "\n";
            // If the last k-mer is in temp, load it fully into k-mer factory
            uint64_t old_temp_slot;
            if (!old_hash_table->hash_table_array[last_position_in_chain].is_in_main_array())
            {
                //std::cout << "was not in main\n";
                // Get the slot in temp
                if (old_hash_table->hash_table_array[chain_position].is_flagged_2())
                {
                    old_temp_slot = hash_table_array[old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot()].get_previous_kmer_slot();
                }
                else
                {
                    old_temp_slot = old_hash_table->hash_table_array[last_position_in_chain].get_previous_kmer_slot();
                }
               
                for (int k = 0; k < kmer_factory->number_of_blocks; k++)
                {
                    //std::cout << "POSSIBLE ERROR SPOT STARTS\n";
                    uint64_t block_to_be_pushed = old_hash_table->temp_array[old_temp_slot*kmer_blocks+k];
                    //std::cout << "GOT OVER IT\n";
                    // If we are pushing the first block do this (possibly not full block)
                    if (k == 0)
                    {
                        //std::cout << "...\n";
                        for (int m = 0; m < kmer_factory->bits_in_last_block; m+=2)
                        {
                            //std::cout << "pushed char: " << (uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2))) << "\n"; 
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (kmer_factory->bits_in_last_block-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                            
                        }
                        //std::cout << "...\n";
                    }
                    // Otherwise we are pushing a full block
                    else 
                    {
                        for (int m = 0; m < 64; m+=2)
                        {
                            kmer_factory->push_new_integer(uint64_t(3) & (block_to_be_pushed >> (64-2)));
                            new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                            block_to_be_pushed = block_to_be_pushed << 2;
                        }
                    }
                }
            } // Now the kmer_factory holds the "first" k-mer completely
            // Alternatively, if the last k-mer is in the main array do this
            else
            {
                //std::cout << "was in main\n";
                at_chain_position = chain_length-1;
                //std::cout << "...\n";
                for (int k = 0; k < kmer_len; k++)
                {
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    at_chain_position-=1;
                }
                //std::cout << "...\n";
                at_chain_position+=1;
            }

            //std::cout << "FULL CHAIN FOUND\n";
            //std::cout << "START ADDING FROM OLD POSITION " << old_positions_chain[at_chain_position] << "\n";

            // Now k-mer is loaded with a k-mer whose position in the old hash table is in chain_positions[at_chain_position]

            // Start checking the k-mers and add them to the new hash table it flag2 is not set
            uint64_t at_chain_position_backup = at_chain_position;
            while(at_chain_position >= 0)
            {
                // Get the position in the old array
                uint64_t old_position = old_positions_chain[at_chain_position];
                //std::cout << "- MIGRATING k-MER FROM OLD POSITION " << old_position << "\n";
                
                // If not flagged2, migrate to new hash table
                if (!old_hash_table->hash_table_array[old_position].is_flagged_2())
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << old_hash_table->hash_table_array[old_position].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- actually transfering position " << old_position << "\n";
                    // Initial new position
                    uint64_t new_position = new_hasher->get_current_hash();
                    uint64_t probing_backup = new_position;
                    // Probe for empty new position
                    while(hash_table_array[new_position].is_occupied())
                    {
                        new_position += probe_hasher->probe_1(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
                        if (new_position >= size)
                            new_position = new_position % size;
                        if (new_position == probing_backup){
                            std::cout << "RESIZING FAILURE, VERY LETHAL\n";
                            exit(1);
                        }
                    }

                    // Insert into the found empty position
                    hash_table_array[new_position].set_data(old_hash_table->hash_table_array[old_position].get_data());
                    hash_table_array[new_position].set_count(old_hash_table->hash_table_array[old_position].get_count());
                    //hash_table_array[new_position].set_previous_kmer_slot(old_position);
                    hash_table_array[new_position].set_previous_kmer_slot(old_hash_table->hash_table_array[old_position].get_previous_kmer_slot());
                    hash_table_array[new_position].unset_is_flagged_1();
                    hash_table_array[new_position].unset_is_flagged_2();
                    hash_table_array[new_position].unset_is_written_in_output();
                    
                    // Remember where the old position was migrated to
                    //where_did_i_go[old_position] = new_position;
                    old_hash_table->hash_table_array[old_position].set_previous_kmer_slot(new_position);
                    // Set flag
                    old_hash_table->hash_table_array[old_position].set_is_flagged_2();
                    migrated_count+=1;

                    if (round == 0)
                        round_one_kmers+=1;
                    else
                        round_two_kmers+=1;

                    // Count new slot
                    // Add
                    // Flags
                }
                else
                {
                    //std::cout << "-- the previous k-mer for this position lies at slot " << hash_table_array[old_hash_table->hash_table_array[old_position].get_previous_kmer_slot()].get_previous_kmer_slot() << "\n";
                    //std::cout << "--- not transfering position " << old_position <<  " because already transferred\n";
                }

                at_chain_position-=1;
                // All is done
                if (at_chain_position == -1)
                {
                    break;
                }
                // Or continue with the next
                else
                {
                    //std::cout << "...\n";
                    //std::cout << "pushed char: " << uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()) << "\n";
                    kmer_factory->push_new_integer(uint64_t(old_hash_table->hash_table_array[old_positions_chain[at_chain_position]].get_character()));
                    new_hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
                    //std::cout << "...\n";
                }
            }
            //old_hash_table->hash_table_array[check_position].set_is_flagged_2();
            check_position+=1;
        }
    }

    // SET THE CORRECT PREVIOUS K-MERS IN THE NEW HASH TABLE
    uint64_t check_position_old;
    uint64_t check_position_old_previous;
    uint64_t check_position_new_previous;
    check_position = 0;
    while (check_position < size)
    {
        // If not occupied, skip
        if (!hash_table_array[check_position].is_occupied())
        {
            check_position += 1;
            continue;
        }
        // If in main array, make the swap
        if (hash_table_array[check_position].is_in_main_array())
        {
            hash_table_array[check_position].set_previous_kmer_slot(old_hash_table->hash_table_array[hash_table_array[check_position].get_previous_kmer_slot()].get_previous_kmer_slot());
        }
        // Otherwise, swap is not necessary
        /*
        hash_table_array[check_position].set_previous_kmer_slot();
        check_position_old = hash_table_array[check_position].get_previous_kmer_slot();
        // If old k-mer was in main
        if (old_hash_table->hash_table_array[check_position_old].is_in_main_array())
        {
            check_position_old_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
            check_position_new_previous = where_did_i_go[check_position_old_previous];
        }
        // And if it was not
        else
        {
            check_position_new_previous = old_hash_table->hash_table_array[check_position_old].get_previous_kmer_slot();
        }
        hash_table_array[check_position].set_previous_kmer_slot(check_position_new_previous);
        }*/
        
        check_position += 1;
    }
    std::cout << "Resizing finished\n";
    std::cout << "k-mers migrated in round one: " << round_one_kmers << "\n";
    std::cout << "k-mers migrated in round two: " << round_two_kmers << "\n";
}


// ==============================================================================================================