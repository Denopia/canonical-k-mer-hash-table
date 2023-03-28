#include "kmer_hash_table.hpp"

// === For CANONICAL pointer hash table =========================================================================================

// s = slots, k = k-mer length, b = 64bit blocks per k-mer
PointerHashTableCanonical::PointerHashTableCanonical(uint64_t s, uint64_t k, uint64_t b)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMer[size];
    bits_per_char = 2;
    inserted_items = 0;
    //inserted_complete_kmers = 0;
    kmer_blocks = b;
    probe_hasher = new ProbeHasher1();
    probing_prime = mathfunctions::next_prime(uint64_t(std::floor(size/13.0)));
    // Secondary array stuff
    max_secondary_slots = 100;
    touched_secondary_slots = 0;
    secondary_slots_in_use = 0;
    max_secondary_slot_in_use = 0;
    smallest_unused_secondary_slot = 0;
    secondary_array = std::vector<uint64_t>(b*max_secondary_slots, uint64_t(0));
    secondary_free_slots = std::vector<uint8_t>(max_secondary_slots, 1);

    max_kmer_reconstruction_chain = 0;
    total_reconstruction_chain = 0;
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

bool PointerHashTableCanonical::slot_is_occupied(uint64_t slot)
{
    return hash_table_array[slot].is_occupied();
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

uint64_t PointerHashTableCanonical::get_number_of_inserted_items_in_main()
{
    return inserted_items-secondary_slots_in_use;
}

uint64_t PointerHashTableCanonical::get_number_of_max_secondary_slots()
{
    return  max_secondary_slots;
}

uint64_t PointerHashTableCanonical::get_number_of_secondary_slots_in_use()
{
    return secondary_slots_in_use;
}

uint64_t PointerHashTableCanonical::get_max_number_of_secondary_slots_in_use()
{
    return max_secondary_slot_in_use;
}



uint64_t PointerHashTableCanonical::process_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // Try to find and increment
    uint64_t kmer_slot = find_and_increment(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    // If that did not work, insert new
    if (kmer_slot == size)
    {
        //std::cout << "** k-mer was not found, inserting it as new **\n";
        //std::cout << "k-mer being inserted is " << kmer_factory->get_canonical_block(0) << "\n";

        kmer_slot = insert_new_kmer(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    }
    else
    {
        //std::cout << "** k-mer was found, increased its count\n";
    }
    return kmer_slot;
}


bool PointerHashTableCanonical::check_for_cycle(uint64_t reconstruction_slot, uint64_t avoid_slot)
{
     if (!hash_table_array[reconstruction_slot].is_occupied())
        return false;
    
    uint64_t position = reconstruction_slot;
    // Leftmost untaken character position
    int L;
    // Rightmost untaken character position
    int R;
    L = 0;
    R = kmer_len - 1;
    // Leftmost chain character wrt the initial L value
    int Lc = 0;
    // Rightmost chain character wrt the initial R value
    int Rc = kmer_len - 1;
    // Process in reverse: handle the current chain k-mer in reverse orientation
    bool pir = false;

    bool perform_secondary_array_check = false;

    while(true)
    {
        if (position == avoid_slot)
            return true;
        //std::cout << "In this iteration L is " << L << " and R is " << R << "\n";
        //std::cout << "pir is " << pir << "\n";
        if (!hash_table_array[position].predecessor_exists())
        {
            //std::cout << "Going to perform secondary array reconstruction\n";
            perform_secondary_array_check = true;
            break;
        }
        // Compare left characters
        if (L == Lc)
        {
            L = L + 1;
            if (L > R)
                break;
        }
        // Compare right characers
        if (R == Rc)
        {
            R = R - 1;
            if (L > R)
                break;
        }
        // Modify chain character positions and process in reverse indicator
        if (hash_table_array[position].canonical_during_insertion_self())
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // T T T
                if (pir)
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
                // T T F
                else 
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
            }
            else
            {
                // T F T
                if (pir)
                {
                    // Extend right, flip 
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
                // T F F
                else 
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
            }
            
        }
        else
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // F T T
                if (pir)
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
                // F T F
                else 
                {
                    // Extend right, flip pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
            }
            else
            {
                // F F T
                if (pir)
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
                // F F F
                else 
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
            }
        }
        position = hash_table_array[position].get_predecessor_slot();
    }
    return false;
}


uint64_t PointerHashTableCanonical::find_and_increment(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    uint64_t kmer_slot = find(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    
    if (kmer_slot != size)
    {
        hash_table_array[kmer_slot].increase_count();

        // Move k-mer from secondary array to main
        if(predecessor_exists && !hash_table_array[kmer_slot].predecessor_exists())
        {   
            //std::cout << "Trying to migrate from secondary\n";
            if (check_for_cycle(predecessor_slot, kmer_slot))
            {
                //std::cout << "Cycle detected, migration aborted\n";
                return kmer_slot;
            }
            secondary_slots_in_use -= 1;
            //std::cout << "No cycle, migrated successfully\n";
            uint64_t slot_in_secondary = hash_table_array[kmer_slot].get_predecessor_slot();
            if (kmer_factory->forward_kmer_is_canonical()){
                // Canonical during insertion
                hash_table_array[kmer_slot].set_canonical_during_insertion_self();
                // Set left and right characters
                hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1)));
                hash_table_array[kmer_slot].set_right_character(kmer_factory->get_forward_char_at_position(kmer_len-1));
            } else {
                // Reverse canonical during insertion
                hash_table_array[kmer_slot].unset_canonical_during_insertion_self();
                // Set characters swapping and reversing them
                hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1)));
                hash_table_array[kmer_slot].set_right_character(kmer_factory->get_backward_char_at_position(kmer_len-1));
            }
            // Set predecessor exists
            hash_table_array[kmer_slot].set_predecessor_exists();
            // Set predecessor slot
            hash_table_array[kmer_slot].set_predecessor_slot(predecessor_slot);
            // Set predecessor canonical orientation
            if (kmer_factory->previous_forward_kmer_was_canonical()){
                hash_table_array[kmer_slot].set_canonical_during_insertion_predecessor();
            } else {
                hash_table_array[kmer_slot].unset_canonical_during_insertion_predecessor();
            }
            // Free slot in secondary array
            secondary_free_slots[slot_in_secondary] = 1;
            for (int o = 0; o < kmer_blocks; o++)
                secondary_array[slot_in_secondary*kmer_blocks+o] = 0;
        }
    }
    return kmer_slot;
}

// Find where a k-mer is in the hash table
uint64_t PointerHashTableCanonical::find(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    //std::cout << "FIND INITIATED\n";
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }
    //std::cout << "WITH SLOT " << initial_position << "\n";

    uint64_t kmer_slot = initial_position;

    uint64_t probe_iteration = 1;

    int quick_result;

    while (hash_table_array[kmer_slot].is_occupied())
    {
        quick_result = 0;
        if (predecessor_exists && hash_table_array[kmer_slot].predecessor_exists())
        {
            quick_result = quick_kmer_slot_check_sus(kmer_factory, kmer_slot, predecessor_slot);
            if (quick_result == 1)
            {
                return kmer_slot;
            }
        }

        if (quick_result != -1)
        {
            if (full_kmer_slot_check(kmer_factory, kmer_slot))
            {
                //if (quick_result == -1)
                //    std::cout << "Quick check fails but full check is success?????????????????????????\n";
                //std::cout << "FULL K-MER SLOT CHECK STARTS\n";
                return kmer_slot;
            }
        }

        // Probe to next position
        if (kmer_factory->forward_kmer_is_canonical())
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        else
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        probe_iteration += 1;
        // Check if we have looped around back to the initial position
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == initial_position){
            std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
            exit(1);
        }
    }
    return size;
}

// return: -1 = sure false, 0 = unsure, 1 = sure match
int PointerHashTableCanonical::quick_kmer_slot_check_sus(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot, uint64_t predecessor_slot)
{
    //return 0;
    uint64_t lchar;
    uint64_t rchar;

    if (kmer_factory->forward_kmer_is_canonical())
    {
        lchar = kmer_factory->get_forward_char_at_position(0);
        rchar = kmer_factory->get_forward_char_at_position(kmer_len-1);
    }
    else
    {
        lchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1));
        rchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(0));
    }
    if (lchar != hash_table_array[kmer_slot].get_left_character())
        return -1;
    if (rchar != hash_table_array[kmer_slot].get_right_character())
        return -1;

    // Then match the pointer
    if (predecessor_slot != hash_table_array[kmer_slot].get_predecessor_slot())
    {
        return 0;
    }

    // Finally match orientations?
    return 1;
}

bool PointerHashTableCanonical::full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
{

    if (!hash_table_array[kmer_slot].is_occupied())
    {
        std::cout << "WTF??\n";
        return false;
    }
        

    uint64_t position = kmer_slot;
    // Leftmost unchecked character position
    int L;
    // Rightmost unchecked character position
    int R;
    L = 0;
    R = kmer_len - 1;
    // Leftmost chain character wrt the initial L value
    int Lc = 0;
    // Rightmost chain character wrt the initial R value
    int Rc = kmer_len - 1;
    // Process in reverse: handle the current chain k-mer in reverse orientation
    bool pir = false;
    //std::cout << "\nSTARTING NEW COMPARISON\n";
    //std::cout << "Initial L is " << L << "\n";
    //std::cout << "Initial R is " << R << "\n";
    //std::cout << "Initial L' is " << Lc << "\n";
    //std::cout << "Initial R' is " << Rc << "\n";

    // If L > R i.e. factory is empty
    if (L > R)
    {
        return true;
    }

    uint64_t lchar;
    uint64_t rchar;

    if (kmer_factory->forward_kmer_is_canonical())
    {
        lchar = kmer_factory->get_forward_char_at_position(L);
        rchar = kmer_factory->get_forward_char_at_position(R);
    }
    else
    {
        lchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-L));
        rchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-R));
    }

    bool perform_secondary_array_check = false;

    while(true)
    {
        //if (pir)
            //std::cout << "  Looking at position " << position << " in reverse\n";
        //else
            //std::cout << "  Looking at position " << position << " in forward\n";
        // Confirm that the k-mer has a predecessor. If not, move to the special algorithm that handles k-mers in secondary array
        
        // If the current k-mer is fully in the secondary array, break out and perform separate check
        if (!hash_table_array[position].predecessor_exists())
        {
            perform_secondary_array_check = true;
            break;
        }

        // Compare left characters
        //std::cout << "Current k-mer factory left char is: " << lchar << "\n";
        //std::cout << "Current k-mer factory right char is: " << rchar << "\n";

        if (L == Lc)
        {
            //std::cout << "Comparing left characters\n";
            // If needs to be handled in reverse
            if (pir)
            {
                //std::cout << "PIR is true\n";
                if (lchar != twobitstringfunctions::reverse_int(hash_table_array[position].get_right_character()))
                {
                    return false;
                }
                else
                {
                    L = L + 1;
                    if (L > R)
                    {
                        return true;
                    }
                }
            }
            // If needs to be handled forward
            else
            {
                //std::cout << "PIR is false\n";
                if (lchar != hash_table_array[position].get_left_character())
                {
                    return false;
                }
                else
                {
                    L = L + 1;
                    if (L > R)
                    {
                        return true;
                    }
                }
            }
            // Update the leftmost unchecked character
            if ((L >= 0) && (L < kmer_len))
            {
                if (kmer_factory->forward_kmer_is_canonical())
                {
                    lchar = kmer_factory->get_forward_char_at_position(L);
                }
                else
                {
                    lchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-L));
                }
            }
        }

        // Compare right characers
        if (R == Rc)
        {
            //std::cout << "Comparing right characters\n";
            // If needs to be handled in reverse
            if (pir)
            {
                //std::cout << "PIR is true\n";
                if (rchar != twobitstringfunctions::reverse_int(hash_table_array[position].get_left_character()))
                {
                    return false;
                }
                else
                {
                    R = R - 1;
                    if (L > R)
                    {
                        return true;
                    }
                }
            }
            // If needs to be handled forward
            else
            {
                //std::cout << "PIR is false\n";
                if (rchar != hash_table_array[position].get_right_character())
                {
                    return false;
                }
                else
                {
                    R = R - 1;
                    if (L > R)
                    {
                        return true;
                    }
                }
            }
            // Update the rightmost unchecked character
            if ((R >= 0) && (R < kmer_len))
            {
                if (kmer_factory->forward_kmer_is_canonical())
                {
                    rchar = kmer_factory->get_forward_char_at_position(R);
                }
                else
                {
                    rchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-R));
                }
            }
        }

        // Modify chain character positions and process in reverse indicator
        if (hash_table_array[position].canonical_during_insertion_self())
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // T T T
                if (pir)
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
                // T T F
                else 
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
            }
            else
            {
                // T F T
                if (pir)
                {
                    // Extend right, flip 
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
                // T F F
                else 
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
            }
            
        }
        else
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // F T T
                if (pir)
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
                // F T F
                else 
                {
                    // Extend right, flip pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
            }
            else
            {
                // F F T
                if (pir)
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
                // F F F
                else 
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
            }
        }
        //std::cout << "Current L is " << L << "\n";
        //std::cout << "Current R is " << R << "\n";
        //std::cout << "Current L' is " << Lc << "\n";
        //std::cout << "Current R' is " << Rc << "\n";
        position = hash_table_array[position].get_predecessor_slot();
    }

    if (perform_secondary_array_check)
    {
        int Ls = L - Lc;
        int Rs = R - Lc;
        //int Rs = Ls + kmer_len - 1; 
        int a;
        int b;
        uint64_t query_char;
        uint64_t array_char;
        uint64_t secondary_array_position = hash_table_array[position].get_predecessor_slot();
        //std::cout << "#########################################################################################\n";
        //std::cout << "Starting secondary array search when\n";
        //std::cout << "L = " << L << " and L' = " << Lc << " and L'' = " << Ls << "\n";
        //std::cout << "R = " << R << " and R' = " << Rc << " and R'' = " << Rs << "\n";
        //std::cout << "Secondary array is " << secondary_array[secondary_array_position] << "\n";
        //std::cout << "#########################################################################################\n";
        if (!pir)
        {
            //std::cout << "doing it forward\n";
            a = L;
            b = Ls;
            while(a <= R)
            {
                if (kmer_factory->forward_kmer_is_canonical())
                    query_char = kmer_factory->get_forward_char_at_position(a);
                else
                    query_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-a));
                //std::cout << "Query char at a = " << a << " is " << query_char << "\n";
                array_char = get_secondary_array_char(secondary_array_position, b);
                //std::cout << "Array char at b = " << b << " is " << array_char << "\n";
                if (query_char != array_char)
                {
                    return false;
                }
                a+=1;
                b+=1;
            }
            return true;
        }
        else
        {  
            //std::cout << "doing it in reverse\n";
            a = L;
            b = kmer_len-Ls-1;
            while(a <= R)
            {
                if (kmer_factory->forward_kmer_is_canonical())
                    query_char = kmer_factory->get_forward_char_at_position(a);
                else
                    query_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-a));
                //std::cout << "Query char at a = " << a << " is " << query_char << "\n";
                array_char = twobitstringfunctions::reverse_int(get_secondary_array_char(secondary_array_position, b));
                //std::cout << "Array char at b = " << b << " is " << array_char << "\n";
                if (query_char != array_char)
                {
                    return false;
                }
                a+=1;
                b-=1;
            }
            return true;
        }
    }
    return false;
}


uint64_t PointerHashTableCanonical::get_secondary_array_char(uint64_t secondary_array_position, int char_position)
{
    //std::cout << "Asking for secondary array position " << secondary_array_position << " character at position " << char_position << "\n";
    if ((char_position < 0) || (char_position > kmer_len - 1))
    {
        std::cout << "Error in checking k-mer in the secondary array\n";
        exit(1);
    }
    uint64_t pos = kmer_len - char_position - 1;
    uint64_t block = pos / 32;
    uint64_t block_offset = kmer_blocks - block - 1;
    uint64_t block_pos = pos % 32;
    uint64_t return_char = secondary_array[secondary_array_position*kmer_blocks + block_offset];
    //std::cout << "Secondary block is " << secondary_array[secondary_array_position*kmer_blocks + block_offset] << "\n";
    return_char = return_char >> (2*block_pos);
    return_char = return_char & uint64_t(3);
    //std::cout << "Returning charatcer " << return_char << "\n"; 
    return return_char;
}

bool PointerHashTableCanonical::full_kmer_slot_check_NO_SECONDARY(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
{
    if (!hash_table_array[kmer_slot].is_occupied())
        return false;

    uint64_t position = kmer_slot;
    // Leftmost unchecked character position
    int L;
    //int L = 0;
    // Rightmost unchecked character position
    int R;
    if ((hash_table_array[position].canonical_during_insertion_self()) || (true))
    {
        L = kmer_len - kmer_factory->get_number_of_stored_characters();
        //int L = 0;
        // Rightmost unchecked character position
        R = kmer_len - 1;
    }
    else
    {
        L = 0;
        //int L = 0;
        // Rightmost unchecked character position
        R = kmer_factory->get_number_of_stored_characters() - 1;
    }
    // Leftmost chain character wrt the initial L value
    int Lc = 0;
    // Rightmost chain character wrt the initial R value
    int Rc = kmer_len - 1;
    // Process in reverse: handle the current chain k-mer in reverse orientation
    bool pir = false;
    //std::cout << "\nSTARTING NEW COMPARISON\n";
    //std::cout << "Initial L is " << L << "\n";
    //std::cout << "Initial R is " << R << "\n";
    //std::cout << "Initial L' is " << Lc << "\n";
    //std::cout << "Initial R' is " << Rc << "\n";

    // If L > R i.e. factory is empty
    if (L > R)
    {
        return true;
    }

    uint64_t lchar;
    uint64_t rchar;

    if (kmer_factory->forward_kmer_is_canonical())
    {
        lchar = kmer_factory->get_forward_char_at_position(L);
        rchar = kmer_factory->get_forward_char_at_position(R);
    }
    else
    {
        //lchar = kmer_factory->get_backward_char_at_position(L);
        //rchar = kmer_factory->get_backward_char_at_position(R);
        lchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-L));
        rchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-R));
    }

    while(true)
    {
        //if (pir)
            //std::cout << "  Looking at position " << position << " in reverse\n";
        //else
            //std::cout << "  Looking at position " << position << " in forward\n";
        // Confirm that the k-mer has a predecessor. If not, move to the special algorithm that handles k-mers in secondary array
        // Not applicable in this version

        // Compare left characters
        //std::cout << "Current k-mer factory left char is: " << lchar << "\n";
        //std::cout << "Current k-mer factory right char is: " << rchar << "\n";

        if (L == Lc)
        {
            //std::cout << "Comparing left characters\n";
            // If needs to be handled in reverse
            if (pir)
            {
                //std::cout << "PIR is true\n";
                if (hash_table_array[position].right_char_is_null())
                {
                    //std::cout << "Encountered shorter k-mer than the one that was being searched 1\n";
                    //return false;
                }
                else
                {
                    if (lchar != twobitstringfunctions::reverse_int(hash_table_array[position].get_right_character()))
                    {
                        return false;
                    }
                    else
                    {
                        L = L + 1;
                        if (L > R)
                        {
                            return true;
                        }
                    }
                }
            }
            // If needs to be handled forward
            else
            {
                //std::cout << "PIR is false\n";
                if (hash_table_array[position].left_char_is_null())
                {
                    //std::cout << "Encountered shorter k-mer than the one that was being searched 2\n";
                    //return false;
                }
                else
                {
                    if (lchar != hash_table_array[position].get_left_character())
                    {
                        return false;
                    }
                    else
                    {
                        L = L + 1;
                        if (L > R)
                        {
                            return true;
                        }
                    }
                }
            }
            // Update the leftmost unchecked character
            if ((L >= 0) && (L < kmer_len))
            {
                if (kmer_factory->forward_kmer_is_canonical())
                {
                    lchar = kmer_factory->get_forward_char_at_position(L);
                }
                else
                {
                    lchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-L));
                    //lchar = kmer_factory->get_backward_char_at_position(L);
                }
            }
        }

        // Compare right characers
        if (R == Rc)
        {
            //std::cout << "Comparing right characters\n";
            // If needs to be handled in reverse
            if (pir)
            {
                //std::cout << "PIR is true\n";
                if (hash_table_array[position].left_char_is_null())
                {
                    //std::cout << "Encountered shorter k-mer than the one that was being searched 3\n";
                    //return false;
                }
                else
                {
                    if (rchar != twobitstringfunctions::reverse_int(hash_table_array[position].get_left_character()))
                    {
                        return false;
                    }
                    else
                    {
                        R = R - 1;
                        if (L > R)
                        {
                            return true;
                        }
                    }
                }
            }
            // If needs to be handled forward
            else
            {
                //std::cout << "PIR is false\n";
                if (hash_table_array[position].right_char_is_null())
                {
                    //std::cout << "Encountered shorter k-mer than the one that was being searched 4\n";
                    //return false;
                }
                else
                {
                    if (rchar != hash_table_array[position].get_right_character())
                    {
                        return false;
                    }
                    else
                    {
                        R = R - 1;
                        if (L > R)
                        {
                            return true;
                        }
                    }
                }
            }
            // Update the rightmost unchecked character
            if ((R >= 0) && (R < kmer_len))
            {
                if (kmer_factory->forward_kmer_is_canonical())
                {
                    rchar = kmer_factory->get_forward_char_at_position(R);
                }
                else
                {
                    rchar = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1-R));
                    //rchar = kmer_factory->get_backward_char_at_position(R);
                }
            }
        }

        // If all characters match, report k-mer found
        
        

        // If no predecessor, report k-mer not found
        if (!hash_table_array[position].predecessor_exists())
        {
            //std::cout << "No predecessor! Return false\n";
            return false;
        }

        // Modify chain character positions and process in reverse indicator
        if (hash_table_array[position].canonical_during_insertion_self())
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // T T T
                if (pir)
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
                // T T F
                else 
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
            }
            else
            {
                // T F T
                if (pir)
                {
                    // Extend right, flip 
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
                // T F F
                else 
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
            }
            
        }
        else
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // F T T
                if (pir)
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
                // F T F
                else 
                {
                    // Extend right, flip pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
            }
            else
            {
                // F F T
                if (pir)
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
                // F F F
                else 
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
            }
        }
        //std::cout << "Current L is " << L << "\n";
        //std::cout << "Current R is " << R << "\n";
        //std::cout << "Current L' is " << Lc << "\n";
        //std::cout << "Current R' is " << Rc << "\n";
        position = hash_table_array[position].get_predecessor_slot();
    }
    return false;
}

uint64_t PointerHashTableCanonical::insert_new_kmer_in_secondary(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher)
{
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }
    uint64_t kmer_slot = initial_position;
    uint64_t probe_iteration = 1;
    // Find the next empty slot
    while (hash_table_array[kmer_slot].is_occupied())
    {
        // Probe to next position
        if (kmer_factory->forward_kmer_is_canonical())
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        else
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        probe_iteration += 1;
        // Check if we have looped around back to the initial position
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == initial_position){
            //std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
            exit(1);
        }
    }
    // empty slot in main found, insert nedessary information
    
    // Set count
    hash_table_array[kmer_slot].set_count(1);
    // Set occupied
    hash_table_array[kmer_slot].set_occupied();
    // No predecessor
    hash_table_array[kmer_slot].unset_predecessor_exists();
    // Set info based on own canonical orientation
    if (kmer_factory->forward_kmer_is_canonical()){
        // Canonical during insertion
        hash_table_array[kmer_slot].set_canonical_during_insertion_self();
    } else {
        // Reverse canonical during insertion
        hash_table_array[kmer_slot].unset_canonical_during_insertion_self();
    }
    // Unset flags just to be safe, should not be necessary
    hash_table_array[kmer_slot].unset_is_flagged_1();
    hash_table_array[kmer_slot].unset_is_flagged_2();

    // then add the full k-mer in the secondary array
    
    // Find the next free slot in secondary
    bool empty_secondary_slot_found = false;
    for (int j = 0; j < max_secondary_slots; j++)
    {
        if (secondary_free_slots[j] == 1)
        {
            //std::cout << "Could be improved...\n";
            smallest_unused_secondary_slot = j;
            empty_secondary_slot_found = true;
            break;
        }
    }
    if (!empty_secondary_slot_found)
    {
        //std::cout << "SECONDARY ARRAY SLOT RESIZING ERROR THAT SHOULD NOT HAPPEN\n";
        exit(1);
    }

    inserted_items += 1;
    secondary_slots_in_use += 1;
    max_secondary_slot_in_use = std::max(max_secondary_slot_in_use, smallest_unused_secondary_slot+1);
    touched_secondary_slots = std::max(touched_secondary_slots, max_secondary_slot_in_use);

    secondary_free_slots[smallest_unused_secondary_slot] = 0;

    hash_table_array[kmer_slot].set_predecessor_slot(smallest_unused_secondary_slot);
    for (int i = 0; i < kmer_factory->number_of_blocks; i++)
    {
        //std::cout << "Secondary block " << i << " is " << kmer_factory->get_canonical_block(i) << "\n";
        secondary_array[smallest_unused_secondary_slot*kmer_factory->number_of_blocks + i] = kmer_factory->get_canonical_block(i);
    }
    //std::cout << "===\n";

    if (smallest_unused_secondary_slot == max_secondary_slots-1)
    {
        //std::cout << "Secondary array resizing needed\n";
        max_secondary_slots *= 2;
        secondary_array.resize(kmer_blocks*max_secondary_slots, 0);
        secondary_free_slots.resize(max_secondary_slots, 1);
    }

    // Return the slot where the k-mer was inserted
    //std::cout << "++++++++++++ Inserting in slot " << kmer_slot << " and setting predecessor slot to NULL\n";
    //std::cout << "++ canonical during insertion self = " << kmer_factory->forward_kmer_is_canonical() << "\n";

    return kmer_slot;
}

uint64_t PointerHashTableCanonical::insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // If no predecessor, the k-mer goes to the secondary array
    if (!predecessor_exists)
    {
        //std::cout << "Insert in secondary\n";
        return insert_new_kmer_in_secondary(kmer_factory, hasher);
    }

    inserted_items += 1;
    //std::cout << "Insert in main\n";

    // Otherwise, add it to the main array

    // Find initial position based on canonical orientation
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }

    uint64_t kmer_slot = initial_position;
    uint64_t probe_iteration = 1;
    
    // Find the next empty slot
    while (hash_table_array[kmer_slot].is_occupied())
    {
        // Probe to next position
        if (kmer_factory->forward_kmer_is_canonical())
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        else
            kmer_slot += probe_hasher->probe_2(probe_iteration);
        probe_iteration += 1;
        // Check if we have looped around back to the initial position
        if (kmer_slot >= size)
            kmer_slot = kmer_slot % size;
        if (kmer_slot == initial_position){
            std::cout << "Hash table was full and the k-mer was not found. Resizing is probably needed...\n";
            exit(1);
        }
    }
    //std::cout << "__________________________________________________________________________________________________________\n";
    //std::cout << "++++++++++++ Inserting in slot " << kmer_slot << " and setting predecessor slot to " << predecessor_slot << "\n";
    //std::cout << "++ canonical during insertion self = " << kmer_factory->forward_kmer_is_canonical() << "\n";
    //std::cout << "++ canonical during insertion predecessor = " << kmer_factory->previous_forward_kmer_was_canonical() << "\n";
    

    // kmer_slot is now empty, insert the new k-mer
    
    // Set count
    hash_table_array[kmer_slot].set_count(1);
    // Set occupied
    hash_table_array[kmer_slot].set_occupied();
    // Set info about predecessor
    if (predecessor_exists){
        // Predecessor exists
        hash_table_array[kmer_slot].set_predecessor_exists();
        // Predecessor slot
        hash_table_array[kmer_slot].set_predecessor_slot(predecessor_slot);
        // Set predecessor canonical orientation
        if (kmer_factory->previous_forward_kmer_was_canonical()){
            hash_table_array[kmer_slot].set_canonical_during_insertion_predecessor();
        } else {
            hash_table_array[kmer_slot].unset_canonical_during_insertion_predecessor();
        }
    } else {
        // No predecessor
        hash_table_array[kmer_slot].unset_predecessor_exists();
    }
    // Set info based on own canonical orientation
    if (kmer_factory->forward_kmer_is_canonical()){
        // Canonical during insertion
        hash_table_array[kmer_slot].set_canonical_during_insertion_self();
        // Set left and right characters
        hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1)));
        hash_table_array[kmer_slot].set_right_character(kmer_factory->get_forward_char_at_position(kmer_len-1));
        // If k-mer was complete, validate both characters
        if (kmer_factory->current_kmer_is_real()){
            hash_table_array[kmer_slot].unset_left_char_is_null();
            hash_table_array[kmer_slot].unset_right_char_is_null();
        } else {
            // If not complete k-mer, left character is garbage/null
            hash_table_array[kmer_slot].set_left_char_is_null();
            hash_table_array[kmer_slot].unset_right_char_is_null();
        }
    } else {
        // Reverse canonical during insertion
        hash_table_array[kmer_slot].unset_canonical_during_insertion_self();
        // Set characters swapping and reversing them
        hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1)));
        hash_table_array[kmer_slot].set_right_character(kmer_factory->get_backward_char_at_position(kmer_len-1));
        // If k-mer was complete, validate both characters
        if (kmer_factory->current_kmer_is_real()){
            hash_table_array[kmer_slot].unset_left_char_is_null();
            hash_table_array[kmer_slot].unset_right_char_is_null();
        } else {
            hash_table_array[kmer_slot].unset_left_char_is_null();
            // If not complete k-mer, right character is garbage/null
            hash_table_array[kmer_slot].set_right_char_is_null();
        }
        //std::cout << "Left char is NULL : " << hash_table_array[kmer_slot].left_char_is_null() << "\n";
        //std::cout << "Left char is " << hash_table_array[kmer_slot].get_left_character() << "\n";
        //std::cout << "Right char is NULL : " << hash_table_array[kmer_slot].right_char_is_null() << "\n";
        //std::cout << "Right char is " << hash_table_array[kmer_slot].get_right_character() << "\n";
        //std::cout << "__________________________________________________________________________________________________________\n";
    }
    // Unset flags just to be safe, should not be necessary
    hash_table_array[kmer_slot].unset_is_flagged_1();
    hash_table_array[kmer_slot].unset_is_flagged_2();

    // Return the slot where the k-mer was inserted
    return kmer_slot;
}

std::string PointerHashTableCanonical::reconstruct_kmer_in_slot(uint64_t slot)
{
    if (!hash_table_array[slot].is_occupied())
        return "UNOCCUPIED";
    
    //std::cout << "\nreconstructing k-mer\n";
    std::vector<uint64_t> kmer_characters(kmer_len, 4);
    
    uint64_t position = slot;
    // Leftmost untaken character position
    int L;
    // Rightmost untaken character position
    int R;
    L = 0;
    R = kmer_len - 1;
    // Leftmost chain character wrt the initial L value
    int Lc = 0;
    // Rightmost chain character wrt the initial R value
    int Rc = kmer_len - 1;
    // Process in reverse: handle the current chain k-mer in reverse orientation
    bool pir = false;

    bool perform_secondary_array_check = false;

    uint64_t looked_kmers = 0;

    while(true)
    {
        //std::cout << "In this iteration L is " << L << " and R is " << R << "\n";
        //std::cout << "pir is " << pir << "\n";
        if (!hash_table_array[position].predecessor_exists())
        {
            //std::cout << "Going to perform secondary array reconstruction\n";
            perform_secondary_array_check = true;
            break;
        }
        // Compare left characters
        if (L == Lc)
        {
            //std::cout << "Adding left char " << L << "\n";
            // If needs to be handled in reverse
            if (pir)
            {
                kmer_characters[L] = twobitstringfunctions::reverse_int(hash_table_array[position].get_right_character());
            }
            // If needs to be handled forward
            else
            {
                kmer_characters[L] = hash_table_array[position].get_left_character();
            }
            //std::cout << "Adding left char in slot " << L << " = " << kmer_characters[L] << "\n";
            L = L + 1;
            if (L > R)
                break;
        }
        // Compare right characers
        if (R == Rc)
        {
            //std::cout << "Adding right char " << R << "\n";
            // If needs to be handled in reverse
            if (pir)
            {
                kmer_characters[R] = twobitstringfunctions::reverse_int(hash_table_array[position].get_left_character());
            }
            // If needs to be handled forward
            else
            {
                kmer_characters[R] = hash_table_array[position].get_right_character();
            } 
            //std::cout << "Adding right char in slot " << R << " = " << kmer_characters[R] << "\n";       
            R = R - 1;
            if (L > R)
                break;
        }
        // Modify chain character positions and process in reverse indicator
        if (hash_table_array[position].canonical_during_insertion_self())
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // T T T
                if (pir)
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
                // T T F
                else 
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
            }
            else
            {
                // T F T
                if (pir)
                {
                    // Extend right, flip 
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
                // T F F
                else 
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
            }
            
        }
        else
        {
            if (hash_table_array[position].canonical_during_insertion_predecessor())
            {
                // F T T
                if (pir)
                {
                    // Extend left, flip pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                    pir = !pir;
                }
                // F T F
                else 
                {
                    // Extend right, flip pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                    pir = !pir;
                }
            }
            else
            {
                // F F T
                if (pir)
                {
                    // Extend left, keep pir
                    Lc = Lc - 1;
                    Rc = Rc - 1;
                }
                // F F F
                else 
                {
                    // Extend right, keep pir
                    Lc = Lc + 1;
                    Rc = Rc + 1;
                }
            }
        }
        position = hash_table_array[position].get_predecessor_slot();
        looked_kmers += 1;
    }

    if (perform_secondary_array_check)
    {
        //std::cout << "Adding characters from secondary array\n";
       
        int Ls = L - Lc;
        //int Rs = R - Lc;
        //int Rs = Ls + kmer_len - 1;
        int a;
        int b;
        uint64_t query_char;
        uint64_t array_char;
        uint64_t secondary_array_position = hash_table_array[position].get_predecessor_slot();
        //std::cout << "Secondary array is " << secondary_array[secondary_array_position] << "\n";
        if (!pir)
        {
            a = L;
            b = Ls;
            while(a <= R)
            {
                kmer_characters[a] = get_secondary_array_char(secondary_array_position, b);
                //std::cout << "Adding left char in slot " << a << " = " << kmer_characters[a] << "\n";
                a+=1;
                b+=1;
            }
        }
        else
        {  
            a = L;
            //b = Rs;
            b = kmer_len-Ls-1;
            while(a <= R)
            {
                kmer_characters[a] = twobitstringfunctions::reverse_int(get_secondary_array_char(secondary_array_position, b));
                //std::cout << "Adding left char in slot " << a << " = " << kmer_characters[a] << "\n";
                a+=1;
                b-=1;
            }
        }
    }
    std::string return_kmer = "";
    for (int i = 0; i < kmer_len; i++)
    {
        return_kmer = return_kmer + twobitstringfunctions::int2char(kmer_characters[i]);
    }
    if (looked_kmers > max_kmer_reconstruction_chain)
    {
        std::cout << "New maximum k-mer reconstruction chain encountered: " << looked_kmers << "\n";
        max_kmer_reconstruction_chain = looked_kmers;
    }

    total_reconstruction_chain += looked_kmers;
    return return_kmer;
}

void PointerHashTableCanonical::write_kmers_on_disk_separately(uint64_t min_abundance, std::string output_path)
{
    std::ofstream output_file(output_path);

    total_reconstruction_chain = 0;
    uint64_t number_of_reconstructions = 0;
    uint64_t check_position = 0;
    while (check_position < size)
    {
        if ((hash_table_array[check_position].is_occupied()) && (hash_table_array[check_position].get_count() >= min_abundance))
        {
            output_file << reconstruct_kmer_in_slot(check_position) << " " << std::to_string(hash_table_array[check_position].get_count()) << "\n";
            number_of_reconstructions += 1;
        }
        check_position+=1;
    }
    output_file.close(); // close file
	output_file.clear(); // clear flags
    std::cout << "Total reconstruction chain: " << total_reconstruction_chain << "\n";
    if (number_of_reconstructions > 0)
        std::cout << "Average reconstruction chain: " << total_reconstruction_chain / number_of_reconstructions << "\n";
}


void PointerHashTableCanonical::write_kmers_on_disk_separately_even_faster(uint64_t min_abundance, std::string output_path)
{
    std::ofstream output_file(output_path);
    // First count/find how many times each k-mer is referenced
    //std::vector<uint8_t> referenced(size, 0);
    for (int i = 0; i < size; i++)
    {
        hash_table_array[i].unset_is_flagged_1();
        if ((hash_table_array[i].is_occupied()) && (hash_table_array[i].predecessor_exists()))
        {
            //referenced[hash_table_array[i].get_predecessor_slot()] += 1;
            hash_table_array[hash_table_array[i].get_predecessor_slot()].set_is_flagged_1();
        }
    }
    uint64_t check_position;
    uint64_t chain_position;
    //bool reversed = false;
    //bool extend_left = true;
    std::string starting_kmer;
    std::string starting_kmer_rev;
    uint64_t iteration_1_kmers = 0;
    uint64_t iteration_0_kmers = 0;

    uint64_t pred_pos;
    for(int iteration = 0; iteration < 2; iteration++)
    {
        check_position = 0;
        while (check_position < size)
        {
            // If first iteration and referenced, skip
            //if ((iteration == 0) && (referenced[check_position] > 0))
            if ((iteration == 0) && (hash_table_array[check_position].is_flagged_1()))
            {
                check_position += 1;
                continue;
            }
            // If not occupied or already written, skip
            if ((!hash_table_array[check_position].is_occupied()) || (hash_table_array[check_position].is_flagged_2()))
            {
                check_position += 1;
                continue;
            }

            // Construct first chain k-mer
            starting_kmer = reconstruct_kmer_in_slot(check_position);
            starting_kmer_rev = purestringfunctions::reverse_string(starting_kmer);
            //purestringfunctions::reverse_this_string(starting_kmer_rev);
            // Write all unwritten k-mers in the chain
            chain_position = check_position;
            // Write first
            if (hash_table_array[chain_position].get_count() >= min_abundance)
            {
                if (iteration == 0)
                    iteration_0_kmers+=1;
                else
                    iteration_1_kmers+=1;
                output_file << starting_kmer << " " << hash_table_array[chain_position].get_count() << "\n";
            }
            hash_table_array[chain_position].set_is_flagged_2();
            // Move to previous
            if (hash_table_array[chain_position].predecessor_exists())
            {
                pred_pos = hash_table_array[chain_position].get_predecessor_slot();
                if (hash_table_array[pred_pos].predecessor_exists())
                {   
                    if (hash_table_array[chain_position].canonical_during_insertion_self())
                    {
                        if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                            // F + F
                            starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer.substr(0,kmer_len-1);
                            starting_kmer_rev =  starting_kmer_rev.substr(1,kmer_len-1) + twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_left_character()));
                        } else {
                            // F + R
                            starting_kmer = twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_right_character())) + starting_kmer.substr(0,kmer_len-1);
                            starting_kmer_rev = starting_kmer_rev.substr(1, kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                            starting_kmer.swap(starting_kmer_rev);
                        }
                    }
                    else
                    {
                        if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                            // R + F
                            starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_left_character()));
                            starting_kmer_rev = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer_rev.substr(0, kmer_len-1);
                            starting_kmer.swap(starting_kmer_rev);
                        } else {
                            // R + R
                            starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                            starting_kmer_rev = twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_right_character())) + starting_kmer_rev.substr(0,kmer_len-1);
                        }
                    }
                }
                chain_position = pred_pos;
            }
            else
            {
                check_position += 1;
                continue;
            }
            while(!hash_table_array[chain_position].is_flagged_2())
            {
                if (!hash_table_array[chain_position].predecessor_exists())
                {
                    if (hash_table_array[chain_position].get_count() >= min_abundance)
                    {
                        if (iteration == 0)
                            iteration_0_kmers+=1;
                        else
                            iteration_1_kmers+=1;
                        output_file << reconstruct_kmer_in_slot(chain_position) << " " << std::to_string(hash_table_array[chain_position].get_count()) << "\n";
                    }
                    hash_table_array[chain_position].set_is_flagged_2();
                    break;
                }
                if (hash_table_array[chain_position].get_count() >= min_abundance)
                {
                    if (iteration == 0)
                        iteration_0_kmers+=1;
                    else
                        iteration_1_kmers+=1;
                    output_file << starting_kmer << " " << hash_table_array[chain_position].get_count() << "\n";
                }
                hash_table_array[chain_position].set_is_flagged_2();
                pred_pos = hash_table_array[chain_position].get_predecessor_slot();
                if (hash_table_array[pred_pos].predecessor_exists())
                {
                    if (hash_table_array[pred_pos].predecessor_exists())
                    {   
                        if (hash_table_array[chain_position].canonical_during_insertion_self())
                        {
                            if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                                // F + F
                                starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer.substr(0,kmer_len-1);
                                starting_kmer_rev =  starting_kmer_rev.substr(1,kmer_len-1) + twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_left_character()));
                            } else {
                                // F + R
                                starting_kmer = twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_right_character())) + starting_kmer.substr(0,kmer_len-1);
                                starting_kmer_rev = starting_kmer_rev.substr(1, kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                                starting_kmer.swap(starting_kmer_rev);
                            }
                        }
                        else
                        {
                            if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                                // R + F
                                starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_left_character()));
                                starting_kmer_rev = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer_rev.substr(0, kmer_len-1);
                                starting_kmer.swap(starting_kmer_rev);
                            } else {
                                // R + R
                                starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                                starting_kmer_rev = twobitstringfunctions::int2char(twobitstringfunctions::reverse_int(hash_table_array[pred_pos].get_right_character())) + starting_kmer_rev.substr(0,kmer_len-1);
                            }
                        }
                    }
                }
                chain_position = pred_pos;
            }
            check_position += 1;
        }
    }
    std::cout << "Iteration 0 k-mers: " << iteration_0_kmers << "\n";
    std::cout << "Iteration 1 k-mers: " << iteration_1_kmers << "\n";
    output_file.close(); // close file
	output_file.clear(); // clear flags
}

void PointerHashTableCanonical::write_kmers_on_disk_separately_faster(uint64_t min_abundance, std::string output_path)
{
    std::ofstream output_file(output_path);
    // First count/find how many times each k-mer is referenced
    std::vector<uint8_t> referenced(size, 0);
    for (int i = 0; i < size; i++)
    {
        if ((hash_table_array[i].is_occupied()) && (hash_table_array[i].predecessor_exists()))
        {
            referenced[hash_table_array[i].get_predecessor_slot()] += 1;
        }
    }
    uint64_t check_position;
    uint64_t chain_position;
    bool reversed = false;
    bool extend_left = true;
    std::string starting_kmer;
    uint64_t pred_pos;
    for(int iteration = 0; iteration < 2; iteration++)
    {
        check_position = 0;
        while (check_position < size)
        {
            // If first iteration and referenced, skip
            if ((iteration == 0) && (referenced[check_position] > 0))
            {
                check_position += 1;
                continue;
            }
            // If not occupied or already written, skip
            if ((!hash_table_array[check_position].is_occupied()) || (hash_table_array[check_position].is_flagged_1()))
            {
                check_position += 1;
                continue;
            }
            // Construct first chain k-mer
            starting_kmer = reconstruct_kmer_in_slot(check_position);
            // Write all unwritten k-mers in the chain
            chain_position = check_position;
            // Write first
            if (hash_table_array[chain_position].get_count() >= min_abundance)
            {
                output_file << starting_kmer << " " << hash_table_array[chain_position].get_count() << "\n";
            }
            hash_table_array[chain_position].set_is_flagged_1();
            // Move to previous
            if (hash_table_array[chain_position].predecessor_exists())
            {
                pred_pos = hash_table_array[chain_position].get_predecessor_slot();
                if (hash_table_array[pred_pos].predecessor_exists())
                {   
                    if (hash_table_array[chain_position].canonical_during_insertion_self())
                    {
                        if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                            starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer.substr(0,kmer_len-1);
                        } else {
                            starting_kmer = starting_kmer.substr(0,kmer_len-1);
                            purestringfunctions::reverse_this_string(starting_kmer);
                            starting_kmer = starting_kmer + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                            reversed = !reversed;
                        }
                    }
                    else
                    {
                        if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                            starting_kmer = starting_kmer.substr(1,kmer_len-1);
                            purestringfunctions::reverse_this_string(starting_kmer);
                            starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer;
                        } else {
                            starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                        }
                    }
                }
                chain_position = pred_pos;
            }
            else
            {
                check_position += 1;
                continue;
            }
            while(!hash_table_array[chain_position].is_flagged_1())
            {
                if (!hash_table_array[chain_position].predecessor_exists())
                {
                    if (hash_table_array[chain_position].get_count() >= min_abundance)
                    {
                        output_file << reconstruct_kmer_in_slot(chain_position) << " " << std::to_string(hash_table_array[chain_position].get_count()) << "\n";
                    }
                    hash_table_array[chain_position].set_is_flagged_1();
                    break;
                }
                if (hash_table_array[chain_position].get_count() >= min_abundance)
                {
                    output_file << starting_kmer << " " << hash_table_array[chain_position].get_count() << "\n";
                }
                hash_table_array[chain_position].set_is_flagged_1();
                pred_pos = hash_table_array[chain_position].get_predecessor_slot();
                if (hash_table_array[pred_pos].predecessor_exists())
                {
                    if (hash_table_array[pred_pos].predecessor_exists())
                    {   
                        if (hash_table_array[chain_position].canonical_during_insertion_self())
                        {
                            if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                                starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer.substr(0,kmer_len-1);
                            } else {
                                starting_kmer = starting_kmer.substr(0,kmer_len-1);
                                purestringfunctions::reverse_this_string(starting_kmer);
                                starting_kmer = starting_kmer + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                                reversed = !reversed;
                            }
                        }
                        else
                        {
                            if (hash_table_array[chain_position].canonical_during_insertion_predecessor()){
                                starting_kmer = starting_kmer.substr(1,kmer_len-1);
                                purestringfunctions::reverse_this_string(starting_kmer);
                                starting_kmer = twobitstringfunctions::int2char(hash_table_array[pred_pos].get_left_character()) + starting_kmer;
                            } else {
                                starting_kmer = starting_kmer.substr(1,kmer_len-1) + twobitstringfunctions::int2char(hash_table_array[pred_pos].get_right_character());
                            }
                        }
                    }
                }
                chain_position = pred_pos;
            }
            check_position += 1;
        }
    }
    output_file.close(); // close file
	output_file.clear(); // clear flags
}


void PointerHashTableCanonical::write_kmers_on_disk_in_blocks(KMerFactoryCanonical2BC* kmer_factory, uint64_t min_abundance, std::string output_path)
{
    /*
    std::ofstream output_file(output_path);

    uint64_t check_position = 0;
    uint64_t chain_position = 0;
    std::vector<uint64_t> current_chain_abundances;
    std::vector<bool> current_chain_canonical_during_insertion;
    //std::vector<bool> current_chain_completeness;
    std::string current_chain_string = "";
    std::string current_chain_string_reverse = "";
    int current_chain_streak = 0;
    int round_one_kmers = 0;
    int round_two_kmers = 0;
    bool starting_kmer_canonical_during_insertion = true;
    bool flip_chars = false;

    std::vector<bool> is_referenced(size,false);
    while(check_position < size)
    {
        if (hash_table_array[check_position].predecessor_exists())
        {
            is_referenced[hash_table_array[check_position].get_predecessor_slot()] = true;
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
        while (check_position < size)
        {
            if (round == 0)
            {
                if ((!hash_table_array[check_position].is_occupied()) || (is_referenced[check_position]))
                {
                    check_position+=1;
                    continue;
                }
            }
            
            current_chain_abundances.clear();
            current_chain_canonical_during_insertion.clear();
            current_chain_string = "";
            current_chain_string_reverse = "";
            current_chain_streak = 0;
            flip_chars = false;
            // is_flagged_1 means it is written in output
            if ((!hash_table_array[check_position].is_flagged_1()) && (hash_table_array[check_position].is_occupied()) && (hash_table_array[check_position].get_count() >= min_abundance))
            {
                chain_position = check_position;
                starting_kmer_canonical_during_insertion = hash_table_array[check_position].canonical_during_insertion_self();

                // Chain together all k-mers that have not been written into output
                while (!hash_table_array[chain_position].is_flagged_1() && (hash_table_array[chain_position].predecessor_exists()))
                {
                    current_chain_streak+=1;
                    current_chain_abundances.push_back(hash_table_array[chain_position].get_count());
                    current_chain_canonical_during_insertion.push_back(hash_table_array[chain_position].canonical_during_insertion_self());
                    if (!flip_chars)
                    {

                    }
                    else
                    {

                    }
                    //current_chain_string = twobitstringfunctions::int2char(uint64_t(hash_table_array[chain_position].get_character())) + current_chain_string;
                    hash_table_array[chain_position].set_is_flagged_1();
                    if(hash_table_array[chain_position].prev_kmer_exists())
                    {
                        if (hash_table_array[chain_position].canonical_during_insertion_self() != hash_table_array[chain_position].canonical_during_insertion_predecessor())
                        {
                            flip_chars = true;
                        }
                        else
                        {
                            flip_chars = false;
                        }
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
    */

}

// ==============================================================================================================