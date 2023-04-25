#include "kmer_hash_table.hpp"
#include "functions_kmer_mod.hpp"

// === For CANONICAL pointer hash table =========================================================================================

// s = slots, k = k-mer length, b = 64bit blocks per k-mer
PointerHashTableCanonicalAF::PointerHashTableCanonicalAF(uint64_t s, uint64_t k, uint64_t b)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMerAtomicFlag[size];
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

    //main_locks = new std::atomic_flag[size];
    secondary_lock.clear();

    //for (int i = 0; i < s; i++)
    //    main_locks[i].clear();
    
    max_kmer_reconstruction_chain = 0;
    total_reconstruction_chain = 0;
}

PointerHashTableCanonicalAF::~PointerHashTableCanonicalAF()
{
    delete[] hash_table_array;
    //delete[] main_locks;
    delete probe_hasher;

}

uint64_t PointerHashTableCanonicalAF::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(true); 
}

bool PointerHashTableCanonicalAF::kmer_in_slot_is_complete(uint64_t slot)
{
    return hash_table_array[slot].is_complete(true); 
}

bool PointerHashTableCanonicalAF::slot_is_occupied(uint64_t slot)
{
    return hash_table_array[slot].is_occupied(true);
}

void PointerHashTableCanonicalAF::resize()
{
    std::cout << "Hash table resizing not implemented yet...\n";
    exit(1);
}

uint64_t PointerHashTableCanonicalAF::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t PointerHashTableCanonicalAF::get_number_of_inserted_items_in_main()
{
    return inserted_items-secondary_slots_in_use;
}

uint64_t PointerHashTableCanonicalAF::get_number_of_max_secondary_slots()
{
    return  max_secondary_slots;
}

uint64_t PointerHashTableCanonicalAF::get_number_of_secondary_slots_in_use()
{
    return secondary_slots_in_use;
}

uint64_t PointerHashTableCanonicalAF::get_max_number_of_secondary_slots_in_use()
{
    return max_secondary_slot_in_use;
}


//
//
// MAIN TO MODIFY (0) : DONE, nothing needs to be locked
// Nothing is locked when we call this
//
uint64_t PointerHashTableCanonicalAF::process_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
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

//
//
// MAIN TO MODIFY (5) : DONE
// Nothing is locked when we go here
//
bool PointerHashTableCanonicalAF::check_for_cycle(uint64_t reconstruction_slot, uint64_t avoid_slot)
{
    //while(main_locks[reconstruction_slot].test_and_set(std::memory_order_acquire));
    bool slot_occupied = hash_table_array[reconstruction_slot].is_occupied(true);
    //main_locks[reconstruction_slot].clear(std::memory_order_release);
    if (!slot_occupied)
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
        
        //while(main_locks[position].test_and_set(std::memory_order_acquire));
        bool slot_has_predecessor = hash_table_array[position].predecessor_exists(true);
        bool slot_canonical_during_insertion = hash_table_array[position].canonical_during_insertion_self(true);
        bool predecessor_canonical_during_insertion = hash_table_array[position].canonical_during_insertion_self(true);
        uint64_t slots_predecessor_slot = hash_table_array[position].get_predecessor_slot(true);
        //main_locks[position].clear(std::memory_order_release);
        
        //if (!hash_table_array[position].predecessor_exists())
        if (!slot_has_predecessor)
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
        //if (hash_table_array[position].canonical_during_insertion_self())
        if (slot_canonical_during_insertion)
        {
            //if (hash_table_array[position].canonical_during_insertion_predecessor())
            if (predecessor_canonical_during_insertion)
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
            //if (hash_table_array[position].canonical_during_insertion_predecessor())
            if (predecessor_canonical_during_insertion)
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
        //position = hash_table_array[position].get_predecessor_slot();
        position = slots_predecessor_slot;
    }
    return false;
}


//
//
// MAIN TO MODIFY (1) : DONE
// Nothing is locked when we go here
//
uint64_t PointerHashTableCanonicalAF::find_and_increment(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // If the returned value is "size", nothing is locked
    // If the returned value is less than size, the "kmer_slot" is also locked
    // ACTUALLY, don't leave anything locked
    uint64_t kmer_slot = find(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    
    if (kmer_slot != size)
    {
        //while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
        hash_table_array[kmer_slot].increase_count(true);
        //main_locks[kmer_slot].clear(std::memory_order_release);

        // Move k-mer from secondary array to main
        if(predecessor_exists && !hash_table_array[kmer_slot].predecessor_exists())
        {   
            //std::cout << "Trying to migrate from secondary\n";
            if (check_for_cycle(predecessor_slot, kmer_slot))
            {
                //std::cout << "Cycle detected, migration aborted\n";
                return kmer_slot;
            }
            // Okay, now we go on to change things
            // Lock the kmer_slot again and also lock the secondary slot
            while(secondary_lock.test_and_set(std::memory_order_acquire));
            //while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
            hash_table_array[kmer_slot].acquire_lock();
            // Confirm that k-mer is still in secondary
            if (hash_table_array[kmer_slot].predecessor_exists(false))
            {
                // The k-mer was moved by some other thread, abort mission
                //main_locks[kmer_slot].clear(std::memory_order_release);
                hash_table_array[kmer_slot].release_lock();
                secondary_lock.clear(std::memory_order_release);
                return kmer_slot;
            }
            // k-mer is still in secondary get its secondary array slot
            uint64_t slot_in_secondary = hash_table_array[kmer_slot].get_predecessor_slot(false);
            // Move on to modify information
            secondary_slots_in_use -= 1;
            if (kmer_factory->forward_kmer_is_canonical()){
                // Canonical during insertion
                hash_table_array[kmer_slot].set_canonical_during_insertion_self(false);
                // Set left and right characters
                hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1)),false);
                hash_table_array[kmer_slot].set_right_character(kmer_factory->get_forward_char_at_position(kmer_len-1),false);
            } else {
                // Reverse canonical during insertion
                hash_table_array[kmer_slot].unset_canonical_during_insertion_self(false);
                // Set characters swapping and reversing them
                hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1)),false);
                hash_table_array[kmer_slot].set_right_character(kmer_factory->get_backward_char_at_position(kmer_len-1),false);
            }
            // Set predecessor exists
            hash_table_array[kmer_slot].set_predecessor_exists(false);
            // Set predecessor slot
            hash_table_array[kmer_slot].set_predecessor_slot(predecessor_slot,false);
            // Set predecessor canonical orientation
            if (kmer_factory->previous_forward_kmer_was_canonical()){
                hash_table_array[kmer_slot].set_canonical_during_insertion_predecessor(false);
            } else {
                hash_table_array[kmer_slot].unset_canonical_during_insertion_predecessor(false);
            }
            // Free slot in secondary array
            secondary_free_slots[slot_in_secondary] = 1;
            for (int o = 0; o < kmer_blocks; o++)
                secondary_array[slot_in_secondary*kmer_blocks+o] = 0;
            // Release locks
            //main_locks[kmer_slot].clear(std::memory_order_release);
            hash_table_array[kmer_slot].release_lock();
            secondary_lock.clear(std::memory_order_release);
        }
    }
    return kmer_slot;
}
//
//
// MAIN TO MODIFY (4) : DONE rolled back
// Nothing is locked when we go here
//
// Find where a k-mer is in the hash table
uint64_t PointerHashTableCanonicalAF::find(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
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

    /*
    while(true)
    {
        while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
        bool slot_is_occupied = hash_table_array[kmer_slot].is_occupied();
        bool slot_predecessor_exists = hash_table_array[kmer_slot].predecessor_exists();
        main_locks[kmer_slot].clear(std::memory_order_release);

        if (!slot_is_occupied)
        {
            return size;
        }

        // If queried k-mer has a predecessor and the k-mer in the current slot also has a predecessor, perform a quick check
        quick_result = 0;
        if (predecessor_exists && slot_predecessor_exists)
        {
            quick_result = quick_kmer_slot_check_sus(kmer_factory, kmer_slot, predecessor_slot);
            if (quick_result == 1)
            {
                return kmer_slot;
            }
        }
        // Perform full check
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
    */
    //*
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
    //*/
    return size;
}

//
//
// MAIN TO MODIFY: DONE
// Nothing is locked when we call this
//
// return: -1 = sure false, 0 = unsure, 1 = sure match
int PointerHashTableCanonicalAF::quick_kmer_slot_check_sus(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot, uint64_t predecessor_slot)
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

    /*
    while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
    uint64_t slot_lchar = hash_table_array[kmer_slot].get_left_character();
    uint64_t slot_rchar = hash_table_array[kmer_slot].get_right_character();
    uint64_t slot_predecessor_slot = hash_table_array[kmer_slot].get_predecessor_slot();
    main_locks[kmer_slot].clear(std::memory_order_release);

    if (lchar != slot_lchar)
        return -1;
    if (rchar != slot_rchar)
        return -1;
    if (predecessor_slot != slot_predecessor_slot)
    {
        return 0;
    }
    */
    //*
    if (lchar != hash_table_array[kmer_slot].get_left_character())
        return -1;
    if (rchar != hash_table_array[kmer_slot].get_right_character())
        return -1;

    // Then match the pointer
    if (predecessor_slot != hash_table_array[kmer_slot].get_predecessor_slot())
    {
        return 0;
    }
    //*/
    // Finally match orientations?
    return 1;
}

//
//
// MAIN TO MODIFY : NOT DONE
// Nothing is locked when we go here
//
bool PointerHashTableCanonicalAF::full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
{

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
        while(secondary_lock.test_and_set(std::memory_order_acquire));
        if (hash_table_array[position].predecessor_exists())
        {
            secondary_lock.clear(std::memory_order_release);
            std::cout << "Some thread moved the last k-mer in chain from secondary to main during another thread's k-mer check. Start check from the beginning\n";
            return full_kmer_slot_check(kmer_factory, kmer_slot);
        }
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
                    secondary_lock.clear(std::memory_order_release);
                    return false;
                }
                a+=1;
                b+=1;
            }
            secondary_lock.clear(std::memory_order_release);
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
                    secondary_lock.clear(std::memory_order_release);
                    return false;
                }
                a+=1;
                b-=1;
            }
            secondary_lock.clear(std::memory_order_release);
            return true;
        }
    }
    return false;
}

// DONE ?
uint64_t PointerHashTableCanonicalAF::get_secondary_array_char(uint64_t secondary_array_position, int char_position)
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

bool PointerHashTableCanonicalAF::full_kmer_slot_check_NO_SECONDARY(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
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

//
//
// MAIN TO CHANGE (3) : DONE !
//
//
uint64_t PointerHashTableCanonicalAF::insert_new_kmer_in_secondary(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher)
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
    while(true)
    {
        //while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
        hash_table_array[kmer_slot].acquire_lock();
        if (!hash_table_array[kmer_slot].is_occupied(false))
        {
            break;
        }
        else
        {
            //main_locks[kmer_slot].clear(std::memory_order_release);
            hash_table_array[kmer_slot].release_lock();
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
    }
    // empty slot in main found and it is also locked, insert necessary information
    
    // Set count
    hash_table_array[kmer_slot].set_count(1,false);
    // Set occupied
    hash_table_array[kmer_slot].set_occupied(false);
    // No predecessor
    hash_table_array[kmer_slot].unset_predecessor_exists(false);
    // Set info based on own canonical orientation
    if (kmer_factory->forward_kmer_is_canonical()){
        // Canonical during insertion
        hash_table_array[kmer_slot].set_canonical_during_insertion_self(false);
    } else {
        // Reverse canonical during insertion
        hash_table_array[kmer_slot].unset_canonical_during_insertion_self(false);
    }
    // Unset flags just to be safe, should not be necessary
    hash_table_array[kmer_slot].unset_is_flagged_1(false);
    hash_table_array[kmer_slot].unset_is_flagged_2(false);

    // then add the full k-mer in the secondary array
    
    // Find the next free slot in secondary
    bool empty_secondary_slot_found = false;
    while(secondary_lock.test_and_set(std::memory_order_acquire));
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

    // Empty secondary slot found and the secondary table is still locked

    inserted_items += 1;
    secondary_slots_in_use += 1;
    max_secondary_slot_in_use = std::max(max_secondary_slot_in_use, smallest_unused_secondary_slot+1);
    touched_secondary_slots = std::max(touched_secondary_slots, max_secondary_slot_in_use);

    secondary_free_slots[smallest_unused_secondary_slot] = 0;

    hash_table_array[kmer_slot].set_predecessor_slot(smallest_unused_secondary_slot,false);
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

    // Free main and secondary locks
    secondary_lock.clear(std::memory_order_release);
    //main_locks[kmer_slot].clear(std::memory_order_release);
    hash_table_array[kmer_slot].release_lock();
    // Return the slot where the k-mer was inserted
    return kmer_slot;
}

//
//
// MAIN TO CHANGE (2) : DONE !
//
//
uint64_t PointerHashTableCanonicalAF::insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
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
    while(true)
    {
        //while(main_locks[kmer_slot].test_and_set(std::memory_order_acquire));
        hash_table_array[kmer_slot].acquire_lock();
        if (!hash_table_array[kmer_slot].is_occupied(false))
        {
            break;
        }
        else
        {
            //main_locks[kmer_slot].clear(std::memory_order_release);
            hash_table_array[kmer_slot].release_lock();
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
    }
    // Now "kmer_slot" contains a free slot and we have it locked
    
    // Set count
    hash_table_array[kmer_slot].set_count(1,false);
    // Set occupied
    hash_table_array[kmer_slot].set_occupied(false);
    // Set info about predecessor
    if (predecessor_exists){
        // Predecessor exists
        hash_table_array[kmer_slot].set_predecessor_exists(false);
        // Predecessor slot
        hash_table_array[kmer_slot].set_predecessor_slot(predecessor_slot,false);
        // Set predecessor canonical orientation
        if (kmer_factory->previous_forward_kmer_was_canonical()){
            hash_table_array[kmer_slot].set_canonical_during_insertion_predecessor(false);
        } else {
            hash_table_array[kmer_slot].unset_canonical_during_insertion_predecessor(false);
        }
    } else {
        // No predecessor
        hash_table_array[kmer_slot].unset_predecessor_exists(false);
    }
    // Set info based on own canonical orientation
    if (kmer_factory->forward_kmer_is_canonical()){
        // Canonical during insertion
        hash_table_array[kmer_slot].set_canonical_during_insertion_self(false);
        // Set left and right characters
        hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1)),false);
        hash_table_array[kmer_slot].set_right_character(kmer_factory->get_forward_char_at_position(kmer_len-1),false);
        // If k-mer was complete, validate both characters
        if (kmer_factory->current_kmer_is_real()){
            hash_table_array[kmer_slot].unset_left_char_is_null(false);
            hash_table_array[kmer_slot].unset_right_char_is_null(false);
        } else {
            // If not complete k-mer, left character is garbage/null
            hash_table_array[kmer_slot].set_left_char_is_null(false);
            hash_table_array[kmer_slot].unset_right_char_is_null(false);
        }
    } else {
        // Reverse canonical during insertion
        hash_table_array[kmer_slot].unset_canonical_during_insertion_self(false);
        // Set characters swapping and reversing them
        hash_table_array[kmer_slot].set_left_character(twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1)),false);
        hash_table_array[kmer_slot].set_right_character(kmer_factory->get_backward_char_at_position(kmer_len-1),false);
        // If k-mer was complete, validate both characters
        if (kmer_factory->current_kmer_is_real()){
            hash_table_array[kmer_slot].unset_left_char_is_null(false);
            hash_table_array[kmer_slot].unset_right_char_is_null(false);
        } else {
            hash_table_array[kmer_slot].unset_left_char_is_null(false);
            // If not complete k-mer, right character is garbage/null
            hash_table_array[kmer_slot].set_right_char_is_null(false);
        }
        //std::cout << "Left char is NULL : " << hash_table_array[kmer_slot].left_char_is_null() << "\n";
        //std::cout << "Left char is " << hash_table_array[kmer_slot].get_left_character() << "\n";
        //std::cout << "Right char is NULL : " << hash_table_array[kmer_slot].right_char_is_null() << "\n";
        //std::cout << "Right char is " << hash_table_array[kmer_slot].get_right_character() << "\n";
        //std::cout << "__________________________________________________________________________________________________________\n";
    }
    // Unset flags just to be safe, should not be necessary
    hash_table_array[kmer_slot].unset_is_flagged_1(false);
    hash_table_array[kmer_slot].unset_is_flagged_2(false);

    // Everything has been set, time to free the lock
    //main_locks[kmer_slot].clear(std::memory_order_release);
    hash_table_array[kmer_slot].release_lock();
    // Return the slot where the k-mer was inserted
    return kmer_slot;
}

std::string PointerHashTableCanonicalAF::reconstruct_kmer_in_slot(uint64_t slot)
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

void PointerHashTableCanonicalAF::write_kmers_on_disk_separately(uint64_t min_abundance, std::string& output_path)
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


void PointerHashTableCanonicalAF::write_kmers_on_disk_separately_even_faster(uint64_t min_abundance, std::string& output_path)
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

void PointerHashTableCanonicalAF::write_kmers_on_disk_separately_faster(uint64_t min_abundance, std::string& output_path)
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


// ==============================================================================================================

BasicAtomicHashTable::BasicAtomicHashTable(uint64_t s, uint64_t k)
{
    size = s;
    kmer_len = k;
    kmer_array = new BasicAtomicKMer[s];
}

BasicAtomicHashTable::~BasicAtomicHashTable()
{
    delete[] kmer_array;
}

// ==============================================================================================================

// ATOMIC FLAG VERSION

BasicAtomicFlagHashTableLong::BasicAtomicFlagHashTableLong(uint64_t s, uint32_t k)
{
    size = s;
    kmer_len = k;
    kmer_bytes = k / 4;
    if (k % 4 != 0)
        kmer_bytes+=1;
    kmer_locks = new std::atomic_flag[s];
    kmer_array = new uint8_t[s*kmer_bytes];
    counts = new uint32_t[s];
    for (int i = 0; i < s; i++)
        kmer_locks[i].clear();
}

BasicAtomicFlagHashTableLong::~BasicAtomicFlagHashTableLong()
{
    delete[] kmer_locks;
    delete[] kmer_array;
    delete[] counts;
}

void BasicAtomicFlagHashTableLong::write_kmers(uint64_t min_abundance, std::string& output_path)
{
    std::ofstream output_file(output_path);

    for(int i = 0; i < size; i++)
    {   
        // Write k-mer in the output file only if its count is at least min_abundance
        if (counts[i] >= min_abundance)
        {
            //std::cout << "Map position " << i << "\n";
            int current_byte_pos = kmer_bytes*i;
            int unwritten_byte_chars = kmer_len % 4;
            if (unwritten_byte_chars == 0)
                unwritten_byte_chars = 4;
            for (int j = 0; j < kmer_len; j++)
            {
                
                //std::cout << "Byte position " << current_byte_pos << " and value is " << std::bitset<8>(kmer_array[current_byte_pos]) << "\n";
                output_file << twobitstringfunctions::int2char_small((kmer_array[current_byte_pos]>>(2*(unwritten_byte_chars-1))) & uint8_t(3));
                unwritten_byte_chars--;
                if (unwritten_byte_chars == 0)
                {
                    unwritten_byte_chars = 4;
                    current_byte_pos++;
                }
                
            }
            output_file << " " << counts[i] << "\n";
        }
    }

    output_file.close();
	output_file.clear();
    
}

// ==============================================================================================================

// ATOMIC VARIABLE VERSION

BasicAtomicVariableHashTableLong::BasicAtomicVariableHashTableLong(uint64_t s, uint32_t k)
{
    size = s;
    kmer_len = k;
    kmer_bytes = k / 4;
    if (k % 4 != 0)
        kmer_bytes+=1;
    kmer_array = new uint8_t[s*kmer_bytes];
    counts = new std::atomic<uint32_t>[s];
}

BasicAtomicVariableHashTableLong::~BasicAtomicVariableHashTableLong()
{
    delete[] kmer_array;
    delete[] counts;
}

void BasicAtomicVariableHashTableLong::write_kmers(uint64_t min_abundance, std::string& output_path)
{
    std::ofstream output_file(output_path);

    for(int i = 0; i < size; i++)
    {   
        // Write k-mer in the output file only if its count is at least min_abundance
        if ((counts[i].load(std::memory_order_acquire)>>1) >= min_abundance)
        {
            //std::cout << "Map position " << i << "\n";
            int current_byte_pos = kmer_bytes*i;
            int unwritten_byte_chars = kmer_len % 4;
            if (unwritten_byte_chars == 0)
                unwritten_byte_chars = 4;
            for (int j = 0; j < kmer_len; j++)
            {
                
                //std::cout << "Byte position " << current_byte_pos << " and value is " << std::bitset<8>(kmer_array[current_byte_pos]) << "\n";
                output_file << twobitstringfunctions::int2char_small((kmer_array[current_byte_pos]>>(2*(unwritten_byte_chars-1))) & uint8_t(3));
                unwritten_byte_chars--;
                if (unwritten_byte_chars == 0)
                {
                    unwritten_byte_chars = 4;
                    current_byte_pos++;
                }
                
            }
            output_file << " " << (counts[i].load(std::memory_order_acquire)>>1) << "\n";
        }
    }

    output_file.close();
	output_file.clear();
    
}

// ==============================================================================================================
// ==============================================================================================================
//
//
//
//
//
//
//                                 ***   POINTER + ATOMIC VARIABLE VERSION   ***
//
//
//
//
//
//
// ==============================================================================================================
// ==============================================================================================================


PointerHashTableCanonicalAV::PointerHashTableCanonicalAV(uint64_t s, uint64_t k, uint64_t b)
{
    size = s;
    kmer_len = k;
    hash_table_array = new OneCharacterAndPointerKMerAtomicVariable[size];
    bits_per_char = 2;
    inserted_items = 0;
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
    secondary_lock.clear();
    
    max_kmer_reconstruction_chain = 0;
    total_reconstruction_chain = 0;
}

PointerHashTableCanonicalAV::~PointerHashTableCanonicalAV()
{
    delete[] hash_table_array;
    delete probe_hasher;

}

uint64_t PointerHashTableCanonicalAV::get_kmer_count_in_slot(uint64_t slot)
{
    return hash_table_array[slot].get_count(); 
}

bool PointerHashTableCanonicalAV::kmer_in_slot_is_complete(uint64_t slot)
{
    return hash_table_array[slot].is_complete(); 
}

bool PointerHashTableCanonicalAV::slot_is_occupied(uint64_t slot)
{
    return hash_table_array[slot].is_occupied();
}

void PointerHashTableCanonicalAV::resize()
{
    std::cout << "Hash table resizing not implemented yet...\n";
    exit(1);
}

uint64_t PointerHashTableCanonicalAV::get_number_of_inserted_items()
{
    return inserted_items;
}

uint64_t PointerHashTableCanonicalAV::get_number_of_inserted_items_in_main()
{
    return inserted_items-secondary_slots_in_use;
}

uint64_t PointerHashTableCanonicalAV::get_number_of_max_secondary_slots()
{
    return  max_secondary_slots;
}

uint64_t PointerHashTableCanonicalAV::get_number_of_secondary_slots_in_use()
{
    return secondary_slots_in_use;
}

uint64_t PointerHashTableCanonicalAV::get_max_number_of_secondary_slots_in_use()
{
    return max_secondary_slot_in_use;
}

// NEW FUNCTION TO PROCESS K-MER
// The previous implementation did not work correctly with multiple threads
uint64_t PointerHashTableCanonicalAV::process_kmer_MT(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // First, find the initial k-mer slot based on canonical orientation
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }
    uint64_t kmer_slot = initial_position;
    uint64_t probe_iteration = 1;
    int quick_result;
    uint64_t return_slot = size;
    bool probe_normally = true;

    // Loop for probing positions for the k-mer

// +++ PROCESSING LOOP STARTS +++

    while (true)
    {
        probe_normally = true;
        bool kmer_was_processed_correctly = false;
        // If slot is free, try to insert there

// +++ FREE SLOT FOUND, TRYING TO INSERT THERE +++

        if (!hash_table_array[kmer_slot].is_occupied())
        {
           
// +++ IF K-MER BEING INSERTED DOES NOT HAVE A PREDECESSOR WE NEED TO FIND ROOM FOR IT IN THE SECONDARY ARRAY FIRST +++
            uint64_t predecessor_for_insertion = predecessor_slot;
            if (!predecessor_exists)
            {
                // First, modify the secondary slot if needed
                // Take lock
                bool empty_secondary_slot_found = false;
                while(secondary_lock.test_and_set(std::memory_order_acquire));
                for (int j = 0; j < max_secondary_slots; j++)
                {
                    if (secondary_free_slots[j] == 1)
                    {
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
                secondary_slots_in_use += 1;
                max_secondary_slot_in_use = std::max(max_secondary_slot_in_use, smallest_unused_secondary_slot+1);
                touched_secondary_slots = std::max(touched_secondary_slots, max_secondary_slot_in_use);
                secondary_free_slots[smallest_unused_secondary_slot] = 0;
                for (int i = 0; i < kmer_factory->number_of_blocks; i++)
                {
                    secondary_array[smallest_unused_secondary_slot*kmer_factory->number_of_blocks + i] = kmer_factory->get_canonical_block(i);
                }
                if (smallest_unused_secondary_slot == max_secondary_slots-1)
                {
                    max_secondary_slots *= 2;
                    secondary_array.resize(kmer_blocks*max_secondary_slots, 0);
                    secondary_free_slots.resize(max_secondary_slots, 1);
                }
                predecessor_for_insertion = smallest_unused_secondary_slot;
                // Free lock
                secondary_lock.clear(std::memory_order_release);
            }

// +++ SET UP VARIABLES FOR INSERTED K-MER +++

            // Set up variables
            uint64_t self_count = 1;
            bool self_occupied = true;
            // predecessor_exists
            // predecessor_for_insertion
            bool pred_canonical_during_insertion = true;
            if (!kmer_factory->previous_forward_kmer_was_canonical())
                pred_canonical_during_insertion = false;
            bool self_canonical_during_insertion;
            uint64_t self_left_char;
            uint64_t self_right_char;
            if (kmer_factory->forward_kmer_is_canonical())
            {
                self_canonical_during_insertion = true;
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_forward_char_at_position(kmer_len-1);
            }
            else
            {
                self_canonical_during_insertion = false;
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_backward_char_at_position(kmer_len-1);
            }
            // Try putting the k-mer in to the hash table slot
            uint64_t expected_data = 0ULL;
            bool insertion_was_success = true;

// +++ TRY TO INSERT THE K-MER INTO THE HASH TABLE

            while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
            //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                expected_data, 
                kmod::modify_for_insertion(expected_data, predecessor_exists, pred_canonical_during_insertion, predecessor_for_insertion, self_canonical_during_insertion, self_left_char, self_right_char), 
                std::memory_order_acq_rel,
                //std::memory_order_release,
                std::memory_order_relaxed))
            {
                // If insertion fails due to the expected value being non-zero, break out of CAS and process k-mer again in the current slot
                if (expected_data != uint64_t(0))
                {
                    insertion_was_success = false;
                    break;
                }
                else
                {
                    std::cout << "Random fail...\n";
                }
            }

            // If inserted successfully
            if (insertion_was_success)
            {
                inserted_items+=1;
                return_slot = kmer_slot;
                kmer_was_processed_correctly = true;

            }
            // If, insert failed
            else
            {
                kmer_was_processed_correctly = false;
                probe_normally = false;
                if (!predecessor_exists)
                {
                    // We must free the secondary slot that was filled at the beginning
                    while(secondary_lock.test_and_set(std::memory_order_acquire));
                    secondary_free_slots[predecessor_for_insertion] = 1;
                    for (int i = 0; i < kmer_factory->number_of_blocks; i++)
                    {
                        secondary_array[predecessor_for_insertion*kmer_factory->number_of_blocks + i] = 0;
                    }
                    secondary_lock.clear(std::memory_order_release);
                }
            }
        }
        // Otherwise, check if the k-mer in the slot is the same we are currently processing

// +++ OCCUPIED SLOT FOUND, CHECKING IF IT IS THE CORRECT ONE +++

        else
        {
            bool this_is_the_correct_slot = false;
            // Quick check result is sotred here, initialize to "unsure"
            quick_result = 0;
            if (predecessor_exists && hash_table_array[kmer_slot].predecessor_exists() && false)
            {
                quick_result = quick_kmer_slot_check_sus(kmer_factory, kmer_slot, predecessor_slot);
                if (quick_result == 1)
                {
                    this_is_the_correct_slot = true;
                }
            }
            // If quick check result was inconclusive, perform a full check
            if (quick_result == 0)
            {
                if (full_kmer_slot_check(kmer_factory, kmer_slot))
                {
                    this_is_the_correct_slot = true;
                }
            }
            // Now, if the current slot is the correct one, simply increase count by one

// +++ CORRECT OCCUPIED SLOT FOUND, INCREASING COUNT BY ONE AND MIGRATING IF NECESSARY +++

            if (this_is_the_correct_slot)
            {

// +++ INCREASE COUNT BY ONE +++

                hash_table_array[kmer_slot].increase_count();

// +++ IF CURRENT K-MER HAS A PREDECESSOR BUT THE ONE IN THE HASH TABLE DOES NOT, MIGRATE POINTER +++

                if(predecessor_exists && !hash_table_array[kmer_slot].predecessor_exists())
                {
                    // First, lock the secondary array
                    while(secondary_lock.test_and_set(std::memory_order_acquire));
                    // If a cycle is found, do not migrate 
                    if (!check_for_cycle(predecessor_slot, kmer_slot))
                    {
                        
                        // Okay, now we move on to change things
                        // Set up needed variables
                        uint64_t slot_in_secondary = hash_table_array[kmer_slot].get_predecessor_slot();
                        uint64_t self_left_char;
                        uint64_t self_right_char;
                        uint64_t expected_data;
                        bool self_forward_canonical;
                        bool pred_forward_canonical;
                        // Determine variable values
                        if (kmer_factory->forward_kmer_is_canonical()){
                            self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1));
                            self_right_char = kmer_factory->get_forward_char_at_position(kmer_len-1);
                            self_forward_canonical = true;
                        }
                        else
                        {
                            self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1));
                            self_right_char = kmer_factory->get_backward_char_at_position(kmer_len-1);
                            self_forward_canonical = false;
                        }
                        if (kmer_factory->previous_forward_kmer_was_canonical())
                        {
                            pred_forward_canonical = true;
                        }
                        else
                        {
                            pred_forward_canonical = false;
                        }
                        // Start looping the modification until migrated
                        bool i_did_the_migration = true;
                        expected_data = hash_table_array[kmer_slot].get_data();
                        while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
                        //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                            expected_data, 
                            kmod::modify_for_migration(expected_data, self_left_char, self_right_char, predecessor_slot, self_forward_canonical, pred_forward_canonical),
                            //std::memory_order_release,
                            std::memory_order_acq_rel,
                            std::memory_order_relaxed))
                        {
                            // If another thread had already migrated the k-mer pointer, do nothing
                            if (hash_table_array[kmer_slot].predecessor_exists())
                            {
                                i_did_the_migration = false;
                                break;
                            }
                        }
                        // Free slot in secondary array
                        if (i_did_the_migration)
                        {
                            secondary_free_slots[slot_in_secondary] = 1;
                            secondary_slots_in_use -= 1;
                            for (int o = 0; o < kmer_blocks; o++)
                                secondary_array[slot_in_secondary*kmer_blocks+o] = 0;
                        }
                    }
                    // Release lock
                    secondary_lock.clear(std::memory_order_release);
                }

// +++ NEW ADDITION : IF BOTH HAVE A PREDECESSOR, SWAP IF CURRENT PREDECESSOR HAS BIGGER COUNT

                else if(predecessor_exists && hash_table_array[kmer_slot].predecessor_exists() && false)
                {
                    if (std::max(2*hash_table_array[hash_table_array[kmer_slot].get_predecessor_slot()].get_count(), uint64_t(20)) < hash_table_array[predecessor_slot].get_count())
                    {
                        if (!check_for_cycle(predecessor_slot, kmer_slot))
                        {
                            //uint64_t loop_counter = 0;
                            bool self_forward_canonical_x;
                            bool pred_forward_canonical_x;
                            // Determine variable values
                            if (kmer_factory->forward_kmer_is_canonical()){
                                self_forward_canonical_x = true;
                            }
                            else
                            {
                                self_forward_canonical_x = false;
                            }
                            if (kmer_factory->previous_forward_kmer_was_canonical())
                            {
                                pred_forward_canonical_x = true;
                            }
                            else
                            {
                                pred_forward_canonical_x = false;
                            }
                            // Swap if the count of the current k-mer predecessor is bigger
                            while (std::max(2*hash_table_array[hash_table_array[kmer_slot].get_predecessor_slot()].get_count(), uint64_t(20)) < hash_table_array[predecessor_slot].get_count())
                            {
                                //loop_counter+=1;
                                //std::cout << loop_counter << "\n";
                                uint64_t expected_data_before_swap = hash_table_array[kmer_slot].get_data();
                                
                                while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
                                //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                                    expected_data_before_swap, 
                                    kmod::modify_predecessor_slot_and_orientations(expected_data_before_swap, predecessor_slot, self_forward_canonical_x, pred_forward_canonical_x),
                                    //std::memory_order_release,
                                    std::memory_order_acq_rel,
                                    std::memory_order_relaxed))
                                {
                                std::cout << "Predecessor swap failed, trying again.\n";
                                break;
                                }
                            }
                        }
                        else
                        {
                            //std::cout << "ois ollu cycli\n";
                        }
                    }
                }
                return_slot = kmer_slot;
                kmer_was_processed_correctly = true;
            }
        }
        // Then, if the k-mer was processed correctly, exit the loop

// +++ IF K-MER WAS PROCESSED FULLY IN THIS ROUND, EXIT THE LOOP +++

        if (kmer_was_processed_correctly)
        {
            break;
        }
        // Otherwise, probe to the next hash table slot

// +++ OTHERWISE, PROBE TO THE NEXT POSITION +++

        else
        {
            if (probe_normally)
            {
                // This probing method is dumb, try something else later once the program works correctly
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
            else
            {
                std::cout << "Need to check the same slot again...\n";
            }
        }
    }

    return return_slot;

}


// DONE
// Main function that is called when we want the hash table to process a new k-mer
// Return the slot where the k-mer resides in after it is processed
//
uint64_t PointerHashTableCanonicalAV::process_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // Try to find and increment
    uint64_t kmer_slot = find_and_increment(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    // If that did not work, insert new
    if (kmer_slot == size)
    {
        kmer_slot = insert_new_kmer(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    }
    return kmer_slot;
}

// DONE
// Function to find the given k-mer in the hash table.
// If found, count is increased by one. Return the slot where the k-mer resides in.
// If not found, return the size of the hash table.
//
uint64_t PointerHashTableCanonicalAV::find_and_increment(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    uint64_t kmer_slot = find(kmer_factory, hasher, predecessor_exists, predecessor_slot);
    
    if (kmer_slot != size)
    {
        hash_table_array[kmer_slot].increase_count();
        
        // Move k-mer pointer from secondary array to main
        if(predecessor_exists && !hash_table_array[kmer_slot].predecessor_exists() && false)
        {
            // First, lock the secondary array
            while(secondary_lock.test_and_set(std::memory_order_acquire));
            // If a cycle is found, do not migrate 
            if (check_for_cycle(predecessor_slot, kmer_slot))
            {
                // Release lock
                secondary_lock.clear(std::memory_order_release);
                return kmer_slot;
            }
            // Okay, now we move on to change things
            // Set up needed variables
            uint64_t slot_in_secondary = hash_table_array[kmer_slot].get_predecessor_slot();
            uint64_t self_left_char;
            uint64_t self_right_char;
            uint64_t expected_data;
            bool self_forward_canonical;
            bool pred_forward_canonical;
            // Determine variable values
            if (kmer_factory->forward_kmer_is_canonical()){
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_forward_char_at_position(kmer_len-1);
                self_forward_canonical = true;
            }
            else
            {
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_backward_char_at_position(kmer_len-1);
                self_forward_canonical = false;
            }
            if (kmer_factory->previous_forward_kmer_was_canonical())
            {
                pred_forward_canonical = true;
            }
            else
            {
                pred_forward_canonical = false;
            }
            // Start looping the modification until migrated
            bool i_did_the_migration = true;
            expected_data = hash_table_array[kmer_slot].get_data();
            while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
            //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                expected_data, 
                kmod::modify_for_migration(expected_data, self_left_char, self_right_char, predecessor_slot, self_forward_canonical, pred_forward_canonical),
                //std::memory_order_release,
                std::memory_order_acq_rel,
                std::memory_order_relaxed))
            {
                // If another thread had already migrated the k-mer pointer, do nothing
                if (hash_table_array[kmer_slot].predecessor_exists())
                {
                    i_did_the_migration = false;
                    break;
                }
            }
            // Free slot in secondary array
            if (i_did_the_migration)
            {
                secondary_free_slots[slot_in_secondary] = 1;
                secondary_slots_in_use -= 1;
                for (int o = 0; o < kmer_blocks; o++)
                    secondary_array[slot_in_secondary*kmer_blocks+o] = 0;
            }
            // Release lock
            secondary_lock.clear(std::memory_order_release);
        }
    }
    return kmer_slot;
}


// DONE
// Find where a k-mer is in the hash table
uint64_t PointerHashTableCanonicalAV::find(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }

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

// DONE
//
// return: -1 = sure false, 0 = unsure, 1 = sure match
int PointerHashTableCanonicalAV::quick_kmer_slot_check_sus(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot, uint64_t predecessor_slot)
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



// DONE
// Makes a full k-mer slot check. If the last k-mer in chain is in secondary and it is modified by another thread, the check has to start from the beginning
// Maybe working?
//
bool PointerHashTableCanonicalAV::full_kmer_slot_check(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
{

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
        if (!hash_table_array[position].predecessor_exists())
        {
            perform_secondary_array_check = true;
            break;
        }

        // Compare left characters
        if (L == Lc)
        {
            // If needs to be handled in reverse
            if (pir)
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
            // If needs to be handled forward
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
            // If needs to be handled in reverse
            if (pir)
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
            // If needs to be handled forward
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

        // Modify chain character positions and process-in-reverse indicator
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
                    // Extend right, flip pir
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
        while(secondary_lock.test_and_set(std::memory_order_acquire));
        if (hash_table_array[position].predecessor_exists())
        {
            secondary_lock.clear(std::memory_order_release);
            std::cout << "Some thread moved the last k-mer in chain from secondary to main during another thread's k-mer check. Start check from the beginning\n";
            return full_kmer_slot_check(kmer_factory, kmer_slot);
        }

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
                    secondary_lock.clear(std::memory_order_release);
                    return false;
                }
                a+=1;
                b+=1;
            }
            secondary_lock.clear(std::memory_order_release);
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
                    secondary_lock.clear(std::memory_order_release);
                    return false;
                }
                a+=1;
                b-=1;
            }
            secondary_lock.clear(std::memory_order_release);
            return true;
        }
    }
    return false;
}


// DONE
// No secondary array lock used, it is assumed that the caller has locked the table beforehand
uint64_t PointerHashTableCanonicalAV::get_secondary_array_char(uint64_t secondary_array_position, int char_position)
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

// DONE
// Checks if a k-mer chain would have a too short cycle after migrating a k-mer to main
//
bool PointerHashTableCanonicalAV::check_for_cycle(uint64_t reconstruction_slot, uint64_t avoid_slot)
{
    //while(main_locks[reconstruction_slot].test_and_set(std::memory_order_acquire));
    bool slot_occupied = hash_table_array[reconstruction_slot].is_occupied();
    //main_locks[reconstruction_slot].clear(std::memory_order_release);
    if (!slot_occupied)
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
        
        //while(main_locks[position].test_and_set(std::memory_order_acquire));
        bool slot_has_predecessor = hash_table_array[position].predecessor_exists();
        bool slot_canonical_during_insertion = hash_table_array[position].canonical_during_insertion_self();
        bool predecessor_canonical_during_insertion = hash_table_array[position].canonical_during_insertion_self();
        uint64_t slots_predecessor_slot = hash_table_array[position].get_predecessor_slot();
        //main_locks[position].clear(std::memory_order_release);
        
        //if (!hash_table_array[position].predecessor_exists())
        if (!slot_has_predecessor)
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
        //if (hash_table_array[position].canonical_during_insertion_self())
        if (slot_canonical_during_insertion)
        {
            //if (hash_table_array[position].canonical_during_insertion_predecessor())
            if (predecessor_canonical_during_insertion)
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
            //if (hash_table_array[position].canonical_during_insertion_predecessor())
            if (predecessor_canonical_during_insertion)
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
        //position = hash_table_array[position].get_predecessor_slot();
        position = slots_predecessor_slot;
    }
    return false;
}




// UNUSED, but kept here as a reminder for the other version
bool PointerHashTableCanonicalAV::full_kmer_slot_check_NO_SECONDARY(KMerFactoryCanonical2BC* kmer_factory, uint64_t kmer_slot)
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



// DONE? might have errors...
// Inserts k-mer into the main array.
//
uint64_t PointerHashTableCanonicalAV::insert_new_kmer(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher, bool predecessor_exists, uint64_t predecessor_slot)
{
    // If no predecessor, the k-mer goes to the secondary array
    if (!predecessor_exists)
    {
        //std::cout << "Insert in secondary\n";
        return insert_new_kmer_in_secondary(kmer_factory, hasher);
    }
    //inserted_items += 1;
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
    bool insertion_was_success = false;
    bool inserted_by_increasing = false;
    // Find the next empty slot and insert. If fails, probe further
    while(true)
    {
        if (!hash_table_array[kmer_slot].is_occupied())
        {
            // Set up variables
            uint64_t self_count = 1;
            bool self_occupied = true;
            bool predecessor_exists = true;
            // predecessor_slot
            bool pred_canonical_during_insertion = true;
            if (!kmer_factory->previous_forward_kmer_was_canonical())
                pred_canonical_during_insertion = false;
            bool self_canonical_during_insertion;
            uint64_t self_left_char;
            uint64_t self_right_char;
            if (kmer_factory->forward_kmer_is_canonical())
            {
                self_canonical_during_insertion = true;
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_backward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_forward_char_at_position(kmer_len-1);
            }
            else
            {
                self_canonical_during_insertion = false;
                self_left_char = twobitstringfunctions::reverse_int(kmer_factory->get_forward_char_at_position(kmer_len-1));
                self_right_char = kmer_factory->get_backward_char_at_position(kmer_len-1);
            }
            // Try putting the k-mer in to the hash table slot
            uint64_t expected_data = 0ULL;
            insertion_was_success = true;
            inserted_by_increasing = false;
            while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
            //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                expected_data, 
                kmod::modify_for_insertion(expected_data, predecessor_exists, pred_canonical_during_insertion, predecessor_slot, self_canonical_during_insertion, self_left_char, self_right_char), 
                std::memory_order_acq_rel,
                //std::memory_order_release,
                std::memory_order_relaxed))
            {
                // If insertion fails and some other k-mer is in the slot instead, check if it is the same
                if (expected_data != uint64_t(0))
                {
                    std::cout << "Some thread inserted a k-mer before this in the same slot.\n";
                    uint64_t megaslots = find_and_increment(kmer_factory, hasher, predecessor_exists, predecessor_slot);
                    // If the slot was taken and the inserted k-mer was wrong, probe further
                    if (megaslots == size)
                    {
                        std::cout << "Wrong k-mer was inserted in the same slot just before.\n";
                        insertion_was_success = false;
                    }
                    else
                    {
                        std::cout << "The same k-mer was inserted just before.\n";
                        kmer_slot = megaslots;
                        inserted_by_increasing = true;
                    }
                    break;
                }
                else
                {
                    std::cout << "Random fail...\n";
                }
            }
        }
        if (insertion_was_success)
        {
            break;
        }
        else
        {
            // Probe to next
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
    }
    if (!inserted_by_increasing)
        inserted_items+=1;
    return kmer_slot;
}


// DONE? would not be surprised if had errors...
// Inserts new k-mer to the secondary array.
//
//
uint64_t PointerHashTableCanonicalAV::insert_new_kmer_in_secondary(KMerFactoryCanonical2BC* kmer_factory, RollingHasherDual* hasher)
{
    // First, modify the secondary slot
    // Take lock
    bool empty_secondary_slot_found = false;
    while(secondary_lock.test_and_set(std::memory_order_acquire));
    for (int j = 0; j < max_secondary_slots; j++)
    {
        if (secondary_free_slots[j] == 1)
        {
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
    secondary_slots_in_use += 1;
    max_secondary_slot_in_use = std::max(max_secondary_slot_in_use, smallest_unused_secondary_slot+1);
    touched_secondary_slots = std::max(touched_secondary_slots, max_secondary_slot_in_use);
    secondary_free_slots[smallest_unused_secondary_slot] = 0;
    for (int i = 0; i < kmer_factory->number_of_blocks; i++)
    {
        secondary_array[smallest_unused_secondary_slot*kmer_factory->number_of_blocks + i] = kmer_factory->get_canonical_block(i);
    }
    if (smallest_unused_secondary_slot == max_secondary_slots-1)
    {
        max_secondary_slots *= 2;
        secondary_array.resize(kmer_blocks*max_secondary_slots, 0);
        secondary_free_slots.resize(max_secondary_slots, 1);
    }
    // Free lock
    secondary_lock.clear(std::memory_order_release);


    // Now do what we do in the main insertion function
    // Find initial position based on canonical orientation
    uint64_t initial_position;
    if (kmer_factory->forward_kmer_is_canonical()){
        initial_position = hasher->get_current_hash_forward();
    } else {
        initial_position = hasher->get_current_hash_backward();
    }

    uint64_t kmer_slot = initial_position;
    uint64_t probe_iteration = 1;
    bool insertion_was_success = false;
    bool inserted_by_increasing = false;
    // Find the next empty slot and insert. If fails, probe further
    while(true)
    {
        if (!hash_table_array[kmer_slot].is_occupied())
        {
            // Set up variables
            uint64_t self_count = 1;
            bool self_occupied = true;
            bool predecessor_exists = false;
            // predecessor_slot = smallest_unused_secondary_slot
            bool pred_canonical_during_insertion = true; // Does not matter, doeas not exist
            bool self_canonical_during_insertion;
            uint64_t self_left_char = 0;// Does not matter
            uint64_t self_right_char = 0;// Does not matter
            if (kmer_factory->forward_kmer_is_canonical())
            {
                self_canonical_during_insertion = true;
            }
            else
            {
                self_canonical_during_insertion = false;
            }
            // Try putting the k-mer in to the hash table slot
            uint64_t expected_data = 0ULL;
            insertion_was_success = true;
            inserted_by_increasing = false;
            while(!hash_table_array[kmer_slot].data.compare_exchange_strong(
            //while(!hash_table_array[kmer_slot].data.compare_exchange_weak(
                expected_data, 
                kmod::modify_for_insertion(expected_data, predecessor_exists, pred_canonical_during_insertion, smallest_unused_secondary_slot, self_canonical_during_insertion, self_left_char, self_right_char), 
                //std::memory_order_release,
                std::memory_order_acq_rel,
                std::memory_order_relaxed))
            {
                // If insertion fails and some other k-mer is in the slot instead, check if it is the same
                if (expected_data != uint64_t(0))
                {
                    std::cout << "Some thread inserted a k-mer before this in the same slot.\n";
                    uint64_t megaslots = find_and_increment(kmer_factory, hasher, false, size);
                    // If the slot was taken and the inserted k-mer was wrong, probe further
                    if (megaslots == size)
                    {
                        std::cout << "Wrong k-mer was inserted in the same slot just before (2).\n";
                        insertion_was_success = false;
                    }
                    else
                    {
                        std::cout << "The same k-mer was inserted just before (2).\n";
                        kmer_slot = megaslots;
                        inserted_by_increasing = true;
                    }
                    break;
                }
                else
                {
                    std::cout << "Random fail...\n";
                }
            }
        }
        if (insertion_was_success)
        {
            break;
        }
        else
        {
            // Probe to next
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
    }
    if (inserted_by_increasing)
    {
        // We must free the secondary slot that was filled at the beginning
        while(secondary_lock.test_and_set(std::memory_order_acquire));
        secondary_free_slots[smallest_unused_secondary_slot] = 1;
        for (int i = 0; i < kmer_factory->number_of_blocks; i++)
        {
            secondary_array[smallest_unused_secondary_slot*kmer_factory->number_of_blocks + i] = 0;
        }
        secondary_lock.clear(std::memory_order_release);
    }
    else
    {
        inserted_items+=1;
    }
    return kmer_slot;
}

// 
// This is done with single thread always, don't touch for now
// 
std::string PointerHashTableCanonicalAV::reconstruct_kmer_in_slot(uint64_t slot)
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
        //std::cout << "New maximum k-mer reconstruction chain encountered: " << looked_kmers << "\n";
        max_kmer_reconstruction_chain = looked_kmers;
    }

    total_reconstruction_chain += looked_kmers;
    return return_kmer;
}

/*
void PointerHashTableCanonicalAV::write_kmers_on_disk_separately(uint64_t min_abundance, std::string& output_path)
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
*/


void PointerHashTableCanonicalAV::write_kmers_on_disk_separately_even_faster(uint64_t min_abundance, std::string& output_path)
{
    std::ofstream output_file(output_path);
    uint64_t kmer_data;
    // First count/find how many times each k-mer is referenced
    //std::vector<uint8_t> referenced(size, 0);
    for (int i = 0; i < size; i++)
    {
        //hash_table_array[i].unset_is_flagged_1();
        kmer_data = hash_table_array[i].get_data();
        //while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_unflagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
        while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_unflagged_2(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
        //while(!hash_table_array[i].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_unflagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
        if ((hash_table_array[i].is_occupied()) && (hash_table_array[i].predecessor_exists()))
        {
            //referenced[hash_table_array[i].get_predecessor_slot()] += 1;
            kmer_data = hash_table_array[i].get_data();
            //while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_1(kmer_data),std::memory_order_release, std::memory_order_relaxed));
            while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_1(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
            //while(!hash_table_array[i].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_flagged_1(kmer_data),std::memory_order_release, std::memory_order_relaxed));
            //hash_table_array[hash_table_array[i].get_predecessor_slot()].set_is_flagged_1();
        }
        else
        {
            kmer_data = hash_table_array[i].get_data();
            //while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_unflagged_1(kmer_data),std::memory_order_release, std::memory_order_relaxed));
            while(!hash_table_array[i].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_unflagged_1(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
            //while(!hash_table_array[i].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_unflagged_1(kmer_data),std::memory_order_release, std::memory_order_relaxed));
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
            kmer_data = hash_table_array[chain_position].get_data();
            //while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
            while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
            //while(!hash_table_array[chain_position].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
            //hash_table_array[chain_position].set_is_flagged_2();
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
                    kmer_data = hash_table_array[chain_position].get_data();
                    //while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
                    while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
                    //while(!hash_table_array[chain_position].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
                    //hash_table_array[chain_position].set_is_flagged_2();
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
                kmer_data = hash_table_array[chain_position].get_data();
                //while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
                while(!hash_table_array[chain_position].data.compare_exchange_strong(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_acq_rel, std::memory_order_relaxed));
                //while(!hash_table_array[chain_position].data.compare_exchange_weak(kmer_data, kmod::modify_to_be_flagged_2(kmer_data),std::memory_order_release, std::memory_order_relaxed));
                //hash_table_array[chain_position].set_is_flagged_2();
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
    //std::cout << "Iteration 0 k-mers: " << iteration_0_kmers << "\n";
    //std::cout << "Iteration 1 k-mers: " << iteration_1_kmers << "\n";
    output_file.close(); // close file
	output_file.clear(); // clear flags
}

/*
void PointerHashTableCanonicalAV::write_kmers_on_disk_separately_faster(uint64_t min_abundance, std::string& output_path)
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
*/
