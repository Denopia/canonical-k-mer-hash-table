#include "kmer_hash_table.hpp"


// === For basic hash table 1 =========================================================================================

BasicHashTable1::BasicHashTable1(uint64_t s, uint64_t k)
{
    size = s;
    kmer_len = k;
    single_kmer_blocks = std::ceil((2.0*kmer_len) / 64.0);
    
    bits_per_char = 2;
    inserted_items = 0;
    inserted_complete_kmers = 0;
    solid_kmers = 0;
    
    hash_table_array_kmers = new uint64_t[size*single_kmer_blocks];
    //hash_table_array_kmers = sdsl::int_vector<>(size*single_kmer_blocks, 0); 

    for (int i = 0; i < size*single_kmer_blocks; i++)
        hash_table_array_kmers[i] = 0;

    hash_table_array_counts = new uint32_t[size];
    for (int i = 0; i < size; i++)
        hash_table_array_counts[i] = 0;
}

BasicHashTable1::~BasicHashTable1()
{
    delete[] hash_table_array_kmers;
    delete[] hash_table_array_counts;
}

uint32_t BasicHashTable1::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array_counts[slot]; 
}

void BasicHashTable1::resize()
{
    std::cout << "Resizing not implemented yet\n";
}

uint64_t BasicHashTable1::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t BasicHashTable1::get_number_of_inserted_complete_kmers()
{
    return inserted_complete_kmers;
}

uint64_t BasicHashTable1::get_number_of_solid_kmers()
{
    return solid_kmers;
}

uint64_t BasicHashTable1::insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    while(hash_table_array_counts[kmer_slot] != 0)
    {
        kmer_slot += 1;
        if (kmer_slot >= size)
            kmer_slot = 0;
        if (kmer_slot == kmer_hash)
        {
            std::cout << "Hash table is full, needs to be resized\n";
            resize();
            exit(1);
        }
    }
    // Set count to 1
    hash_table_array_counts[kmer_slot] = 1;
    // Insert the k-mer
    for (int i = 0; i < single_kmer_blocks; i++)
    {
        hash_table_array_kmers[(single_kmer_blocks*kmer_slot)+i] = kmer_factory->blocks[i];
    }

    inserted_items += 1;
    inserted_complete_kmers += 1;

    return kmer_slot;
}

uint64_t BasicHashTable1::find(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    bool looks_like_the_kmer_is_here;
    bool kmer_was_found = false;
    // Look through all consecutive used slots to see if the seeked k-mer exists
    while(hash_table_array_counts[kmer_slot] != 0)
    {
        looks_like_the_kmer_is_here = true;
        // Check if the blocks match between the seeked k-mer in the k-mer factory and the cirrent hash table slot
        for (uint64_t i = 0; i < single_kmer_blocks; i++)
        {
            if (hash_table_array_kmers[single_kmer_blocks*kmer_slot+i] != kmer_factory->blocks[i])
            {
                looks_like_the_kmer_is_here = false;
                break;
            }
        }
        if (looks_like_the_kmer_is_here)
        {
            return kmer_slot;
        }
        else
        {
            kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = 0;
            if (kmer_slot == kmer_hash)
            {
                std::cout << "Hash table fully checked...\n";
                return size;
            }
        }
    }
    return size;
}

uint64_t BasicHashTable1::find_and_increment(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    bool looks_like_the_kmer_is_here;
    bool kmer_was_found = false;
    // Look through all consecutive used slots to see if the seeked k-mer exists
    while(hash_table_array_counts[kmer_slot] != 0)
    {
        looks_like_the_kmer_is_here = true;
        // Check if the blocks match between the seeked k-mer in the k-mer factory and the cirrent hash table slot
        for (uint64_t i = 0; i < single_kmer_blocks; i++)
        {
            if (hash_table_array_kmers[single_kmer_blocks*kmer_slot+i] != kmer_factory->blocks[i])
            {
                looks_like_the_kmer_is_here = false;
                break;
            }
        }
        if (looks_like_the_kmer_is_here)
        {
            hash_table_array_counts[kmer_slot] += 1;
            if (hash_table_array_counts[kmer_slot] == 2)
            {
                solid_kmers += 1;
            }
            kmer_was_found = true;
            break;
        }
        else
        {
            kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = 0;
            if (kmer_slot == kmer_hash)
            {
                std::cout << "Hash table is full, needs to be resized\n";
                resize();
                exit(1);
            }
        }
    }
    if (kmer_was_found)
    {
        return kmer_slot;
    }
    else
    {
        return size;
    }
}

uint64_t BasicHashTable1::find_and_increment_andifnot_insert(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    bool looks_like_the_kmer_is_here;
    bool kmer_was_found = false;
    // Look through all consecutive used slots to see if the seeked k-mer exists
    while(hash_table_array_counts[kmer_slot] != 0)
    {
        looks_like_the_kmer_is_here = true;
        // Check if the blocks match between the seeked k-mer in the k-mer factory and the cirrent hash table slot
        for (uint64_t i = 0; i < single_kmer_blocks; i++)
        {
            if (hash_table_array_kmers[single_kmer_blocks*kmer_slot+i] != kmer_factory->blocks[i])
            {
                looks_like_the_kmer_is_here = false;
                break;
            }
        }
        if (looks_like_the_kmer_is_here)
        {
            hash_table_array_counts[kmer_slot] += 1;
            if (hash_table_array_counts[kmer_slot] == 2)
            {
                solid_kmers += 1;
            }
            kmer_was_found = true;
            break;
        }
        else
        {
            kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = 0;
            if (kmer_slot == kmer_hash)
            {
                std::cout << "Hash table is full, needs to be resized\n";
                resize();
                exit(1);
            }
        }
    }
    if (kmer_was_found)
    {
        return kmer_slot;
    }
    else
    {
        return insert_new_kmer(kmer_factory, kmer_slot);
    }
}
// ==============================================================================================================



// === For pointer hash table 1 =========================================================================================


PointerHashTable1::PointerHashTable1(uint64_t s, uint64_t k)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMer[size];
    bits_per_char = 2;
    inserted_items = 0;
    inserted_complete_kmers = 0;
    solid_kmers = 0;
    probe_hasher = new ProbeHasher1();
    probing_prime = mathfunctions::next_prime(uint64_t(std::floor(size/13.0)));

}

PointerHashTable1::~PointerHashTable1()
{
    delete[] hash_table_array;
    delete probe_hasher;
}

uint32_t PointerHashTable1::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(); 
}

void PointerHashTable1::resize()
{
    std::cout << "Hash table resizing not implemented yet...\n";
    exit(1);
}

uint64_t PointerHashTable1::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t PointerHashTable1::get_number_of_inserted_complete_kmers()
{
    return inserted_complete_kmers;
}

uint64_t PointerHashTable1::get_number_of_solid_kmers()
{
    return solid_kmers;
}

// Check if the currently loaded k-mer is in the given slot
bool PointerHashTable1::kmer_slot_check(KMerFactory2BC* kmer_factory, uint64_t slot)
{
    // If slot is empty, k-mer is not here
    if (!hash_table_array[slot].is_occupied())
    {
        //std::cout << "Was not occupied\n";
        return false;
    }
        
    // If searched k-mer is complete but the k-mer in hash table is not, k-mer is not here
    if ((kmer_factory->get_number_of_stored_characters() >= kmer_len) && (!hash_table_array[slot].is_complete()))
    {
        //std::cout << "Was not complete k-mer\n";
        return false;
    }
        
    // If we pass the initial checks, start comparing characters one by one
    uint64_t dbi = kmer_factory->number_of_blocks-1; // which data block we are looking at
    uint64_t db = kmer_factory->blocks[dbi]; // store current data block content here
    uint64_t cbc = 0; // number of checked data block characters
    uint64_t mc = 0; // matched characters
    for (uint64_t i = 0; i < kmer_factory->get_number_of_stored_characters(); i++)
    {
        // If the current hash table slot is unoccupied the k-mer does not exist
        if (!hash_table_array[slot].is_occupied())
            return false;
        // Check that the characters match
        if (uint64_t(hash_table_array[slot].get_character()) != (db&uint64_t(3))){
            return false;
        // Update to look at next characters if the current ones match
        } else {
            mc+=1;
            if (mc == kmer_factory->get_number_of_stored_characters())
                return true;
            if (hash_table_array[slot].prev_kmer_exists()){
                slot = hash_table_array[slot].get_previous_kmer_slot();
            } else {
                return false;
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
    // If the characters have matched until now, the k-mer exists in the hash table
    return true;
}

// Find where a k-mer is in the hash table
uint64_t PointerHashTable1::find(KMerFactory2BC* kmer_factory, uint64_t kmer_hash)
{
    uint64_t kmer_slot = kmer_hash;
    while (hash_table_array[kmer_slot].is_occupied())
    {
        if (kmer_slot_check(kmer_factory, kmer_slot)){
            return kmer_slot;
        } else {
            kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
            //kmer_slot += 1;
            if (kmer_slot >= size)
                kmer_slot = kmer_slot % size;
            if (kmer_slot == kmer_hash){
                std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
                resize();
                exit(1);
            }
        }        
    }
    return size;
}

// Find where a k-mer is in the hash table
uint64_t PointerHashTable1::find_using_prev(uint32_t new_chracter, uint64_t kmer_slot, uint64_t previous_kmer_slot)
{
    if (hash_table_array[kmer_slot].is_complete())
    {
        if ((previous_kmer_slot != size) & hash_table_array[previous_kmer_slot].is_occupied())
        {
            if (uint64_t(hash_table_array[kmer_slot].get_character()) == new_chracter)
            {
                if(uint64_t(hash_table_array[kmer_slot].get_previous_kmer_slot()) == previous_kmer_slot)
                {
                    return kmer_slot;
                }
            }
        }
    }
    return size;
}



// Insert a new k-me to the hash table, we assume the k-mer is not already in the hash table
uint64_t PointerHashTable1::insert_new_kmer(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t prev_kmer_slot)
{
    // First find the next empty slot
    uint64_t kmer_slot = kmer_hash;
    
    while(hash_table_array[kmer_slot].is_occupied())
    {
        kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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
        hash_table_array[kmer_slot].unset_prev_kmer_exists();
    }
    // Set count to one
    hash_table_array[kmer_slot].set_count(1);
    // Set last character
    hash_table_array[kmer_slot].set_character(kmer_factory->get_newest_character());
    // Set occupied
    hash_table_array[kmer_slot].set_occupied();
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
    return kmer_slot;
}

uint64_t PointerHashTable1::find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = find_using_prev(kmer_factory->get_newest_character(), kmer_hash, previous_kmer_slot);
    if (kmer_slot < size)
    {
        hash_table_array[kmer_slot].increase_count();
        //if (hash_table_array[hash_table_array[kmer_slot].get_previous_kmer_slot()].get_count() < hash_table_array[previous_kmer_slot].get_count())
        //{
        //    hash_table_array[kmer_slot].set_previous_kmer_slot(previous_kmer_slot);
        //}
    }
    else
    {
        kmer_slot = find(kmer_factory, kmer_hash);
        if (kmer_slot < size)
        {
            hash_table_array[kmer_slot].increase_count();
        }
    } 
    return kmer_slot;
}


uint64_t PointerHashTable1::find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = find_and_increment(kmer_factory, kmer_hash, previous_kmer_slot);
    if (kmer_slot == size)
    {
        kmer_slot = insert_new_kmer(kmer_factory, kmer_hash, previous_kmer_slot);
    }
    return kmer_slot;
}



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
    probe_hasher = new ProbeHasher1();
    probing_prime = mathfunctions::next_prime(uint64_t(std::floor(size/13.0)));
    max_temp_slots = std::ceil(size/200);
    touched_temp_slots = 0;
    temp_slots_in_use = 0;
    max_temp_slot_in_use = 0;
    smallest_unused_temp_slot = 0;
    temp_array = std::vector<uint64_t>(b*max_temp_slots,uint64_t(0));
    temp_free_slots = sdsl::bit_vector(max_temp_slots, 1);
}

PointerHashTable2::~PointerHashTable2()
{
    delete[] hash_table_array;
    delete probe_hasher;
}

uint32_t PointerHashTable2::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(); 
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
    
    // If slot is occupied a k-mer in temp array, return false // DOES NOT WORK CORRECTLY IF TWO SEQUENCES START WITH THE EXACT SAME K-MER, FIX LATER
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
    uint32_t temp_array_slot = 0;
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
uint64_t PointerHashTable2::find_and_modify(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = kmer_hash;
    while (hash_table_array[kmer_slot].is_occupied())
    {
        if (previous_kmer_exists)
        {
            if (uint64_t(hash_table_array[kmer_slot].get_character()) == uint64_t(kmer_factory->get_newest_character()))
            {
                if (previous_kmer_slot == hash_table_array[kmer_slot].get_previous_kmer_slot())
                {
                    return kmer_slot;
                }
            }
        }
        int check_result = kmer_slot_check(kmer_factory, kmer_slot);
        if (check_result == 1)
        {
            return kmer_slot;
        } 
        else if ((check_result == -1))
        {
            //std::cout << "Probing to next slot\n";
            kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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
            if (previous_kmer_exists)
            {
                uint32_t temp_array_slot_in_question = hash_table_array[kmer_slot].get_previous_kmer_slot();
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
                hash_table_array[kmer_slot].increase_count();
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

                return kmer_slot;
            }
            else
            {
                hash_table_array[kmer_slot].increase_count();
                return kmer_slot;
            }
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
            kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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
        kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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
        kmer_slot += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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
    for (int j = 0; j < max_temp_slots; j++)
    {
        if (temp_free_slots[j] == 1)
        {
            smallest_unused_temp_slot = j;
            break;
        }
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
        exit(1);
        max_temp_slots *= 1.5;
        temp_array.resize(max_temp_slots, 0);
        //temp_occupied_vector.resize(max_temp_slots);

    }
    //std::cout << "INSERTED IN TEMP USING SLOT " << kmer_slot << " IN MAIN AND SLOT " << used_temp_slots-1 << " IN TEMP\n";
    return kmer_slot;
}

uint64_t PointerHashTable2::find_and_increment(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    uint64_t kmer_slot = find_and_modify(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot);
    if (kmer_slot < size)
        hash_table_array[kmer_slot].increase_count();
    return kmer_slot;
}


uint64_t PointerHashTable2::find_and_increment_andifnot_insert(KMerFactory2BC* kmer_factory, uint64_t kmer_hash, bool previous_kmer_exists, uint64_t previous_kmer_slot)
{
    // Try to find and increment
    uint64_t kmer_slot = find_and_increment(kmer_factory, kmer_hash, previous_kmer_exists, previous_kmer_slot);
    // If that did not work, insert new
    if (kmer_slot == size)
    {
        // If previous k-mer exists, insert normally
        if (previous_kmer_exists)
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

void PointerHashTable2::clean_up_temp_array(KMerFactory2BC* kmer_factory, RollingHasher1* hasher)
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
        
        std::cout << "K-MERS BROUGHT BACK HOME " << brought_back_home << "\n";
        //std::cout << "SLOT " << i << " IN MAIN WAS STORED IN TEMP, BRINGING IT BACK HOME\n";
        brought_back_home+=1;
        temp_slots_in_use -= 1;
        kmer_factory->reset();
        hasher->reset();

        uint64_t handled_characters = 0;
        uint64_t previous_kmer_slot = 0;
        uint64_t slot_for_me = 0;
        uint32_t current_slot_in_temp = hash_table_array[i].get_previous_kmer_slot();

        while (handled_characters < kmer_len-1)
        {
            kmer_factory->push_new_character(get_oldest_char_in_temp_slot_and_shift_left(current_slot_in_temp, kmer_factory));
            hasher->update_rolling_hash(kmer_factory->get_newest_character(), kmer_factory->get_pushed_off_character());
            slot_for_me = hasher->get_current_hash();
            while(hash_table_array[slot_for_me].is_occupied())
            {
                slot_for_me += probe_hasher->probe(kmer_factory->blocks[kmer_factory->number_of_blocks-1], probing_prime);
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

char PointerHashTable2::get_oldest_char_in_temp_slot_and_shift_left(uint32_t slot, KMerFactory2BC* kmer_factory)
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


// ==============================================================================================================