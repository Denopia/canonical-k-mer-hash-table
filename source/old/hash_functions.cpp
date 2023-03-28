#include "hash_functions.hpp"


RollingHasherDual::RollingHasherDual(uint64_t q, uint64_t m)
{
    // q = mod value
    // m = max power
    // d = alphabet size (for multiplication)
    // di = modular multiplicative inverse of d (for division)
    this->q = q;
    this->m = m;
    this->d = 4;
    this->di = modular_multiplicative_inverse(d, q);
    bpc = 2;
    character_mask = uint64_t(3);
    current_hash_forward = 0;
    current_hash_backward = 0;
    hashed_count = 0;
    h = 1;
    for (int i = 0; i < m-1; i++)
        h = (h*d)%q;
}

void RollingHasherDual::update_rolling_hash_in(uint64_t in)
{
    current_hash_forward = (d*current_hash_forward + in) % q;
    get_current_hash_backward = (di*current_hash_backward + reverse_int(in)*h) % q;
    hashed_count += 1;
}

// Should be fine now
void RollingHasherDual::update_rolling_hash_in_and_out(uint64_t in, uint64_t out)
{
    // First update forward
    uint64_t baseline = (d*current_hash_forward + in) % q;
    uint64_t to_be_removed = (d*h*out) % q;
    if (to_be_removed > baseline)
    {
        current_hash_forward = ((baseline + q) - to_be_removed);
    }
    else
    {
        current_hash_forward = baseline - to_be_removed;
    }
    // Next update backward
    baseline = current_hash_backward;
    to_be_removed = reverse_int(out);
    if (to_be_removed > baseline)
    {
        baseline = ((baseline + q) - to_be_removed);
    }
    else
    {
        baseline = baseline - to_be_removed;
    }
    // Simulate division with inverse multiplication
    baseline = (baseline * di) % q;
    current_hash_backward = ((reverse_int(in)*h) + baseline) % q;

    hashed_count += 1;
}

void RollingHasherDual::update_rolling_hash(uint64_t in, uint64_t out)
{
    // Update number of characters in the hash value
    hashed_count = std::min(m, hashed_count+1);
    // Based on the outgoing character and hash count, update hash value accordingly
    if ((out == 0) || (hashed_count < m))
    {
        update_rolling_hash_in(in);
    }
    else
    {
        update_rolling_hash_in_and_out(in, out);
    }
}

uint64_t RollingHasherDual::get_current_hash_forward()
{
    return current_hash_forward;
}

uint64_t RollingHasherDual::get_current_hash_backward()
{
    return current_hash_backward;
}

void RollingHasherDual::reset()
{
    current_hash_forward = 0;
    current_hash_backward = 0;
    hashed_count = 0;
}

// Only factory loads remain to be implemented plus check thatthe factories are looking good with the new idea
void RollingHasherDual::load_full_factory_canonical(KMerFactoryCanonical2BC * kmer_factory)
{
    int current_block_position = 0;
    uint64_t current_block_forward = kmer_factory->blocks_forward[current_block_position];
    int bits_already_read_from_current_block = 64-kmer_factory->bits_in_last_block;
    current_block_forward = current_block_forward << (bits_already_read_from_current_block);

    for (int i = 0; i < kmer_factory->get_number_of_stored_characters(); i++)
    {
        update_rolling_hash_in(((current_block >> (64 - 2)) & uint64_t(3)));
        current_block_forward = current_block_forward << 2;
        bits_already_read_from_current_block += 2;
        if (bits_already_read_from_current_block == 64)
        {
            current_block_position += 1;
            current_block_forward = kmer_factory->blocks_forward[current_block_position];
            bits_already_read_from_current_block = 0;
        }
    }
}

 
RollingHasher1::RollingHasher1(uint64_t q, uint64_t m)
{
    // q = mod value
    // m = max power
    this->q = q;
    this->m = m;
    d = 4;
    h = 1;
    for (int i = 0; i < m-1; i++)
        h = (h*d)%q;
    current_hash = 0;
    hashed_count = 0;
    bpc = 2;
    character_mask = uint64_t(3);
    // Calculate big power
    //AP = A;
    //for (int i = 0; i < A-1; i++)
    //    AP = AP * A;
    //AP = AP % P;
    //P = 331;
    // a bigger one could be 11113 
}

uint64_t RollingHasher1::update_rolling_hash_in(uint64_t in)
{
    current_hash = (d*current_hash + in) % q;
    hashed_count += 1;
    return current_hash;
}

uint64_t RollingHasher1::update_rolling_hash_in_and_out(uint64_t in, uint64_t out)
{
    uint64_t baseline = d*current_hash + in;
    uint64_t to_be_removed = d*h*out;
    if (baseline == to_be_removed){
        //std::cout << "\n\n GOT TO SPECIAL CASE: in and out equal \n\n";
        current_hash = 0;
    } else if (baseline > to_be_removed){
        //std::cout << "HERE IS STANDARD CASE XOX\n";
        current_hash = (baseline - to_be_removed) % q;
    } else {
        //std::cout << "HERE IS SPECIAL CASE 2\n";
        current_hash = q - ((to_be_removed - baseline) % q);
        if (current_hash == q)
        {
            //std::cout << "HERE IS EXTRA SPECIAL CASE\n";
            current_hash = 0;
        }       
    }
    //current_hash = (d*(current_hash - h*out) + in) % q;
    hashed_count += 1;
    return current_hash;
}

void RollingHasher1::load_full_factory(KMerFactory2BC* kmer_factory)
{
    //std::cout << "LOADING FROM FULL FACTORY\n";
    int current_block_position = 0;
    uint64_t current_block = kmer_factory->blocks[current_block_position];
    //std::cout << "ORIGINAL BLOCK IS " << current_block << "\n";
    int bits_already_read_from_current_block = 64-kmer_factory->bits_in_last_block;
    //std::cout << "ALREADY LOADED BITS: " << bits_already_read_from_current_block << "\n";
    current_block = current_block << (bits_already_read_from_current_block);

    for (int i = 0; i < kmer_factory->get_number_of_stored_characters(); i++)
    {
        //std::cout << "CURRENT BLOCK IS " << current_block << "\n";
        //std::cout << "UPDATING WITH " << (((current_block >> (64 - 2)) & uint64_t(3))) << "\n";
        update_rolling_hash_in(((current_block >> (64 - 2)) & uint64_t(3)));
        current_block = current_block << 2;
        bits_already_read_from_current_block += 2;
        hashed_count += 1;
        if (bits_already_read_from_current_block == 64)
        {
            current_block_position += 1;
            current_block = kmer_factory->blocks[current_block_position];
            bits_already_read_from_current_block = 0;
        }
    }
}

void RollingHasher1::load_full_factory_canonical(KMerFactoryCanonical2BC * kmer_factory)
{
    //std::cout << "LOADING FROM FULL FACTORY\n";
    int current_block_position = 0;
    uint64_t current_block = kmer_factory->blocks[current_block_position];
    //std::cout << "ORIGINAL BLOCK IS " << current_block << "\n";
    int bits_already_read_from_current_block = 64-kmer_factory->bits_in_last_block;
    //std::cout << "ALREADY LOADED BITS: " << bits_already_read_from_current_block << "\n";
    current_block = current_block << (bits_already_read_from_current_block);

    for (int i = 0; i < kmer_factory->get_number_of_stored_characters(); i++)
    {
        //std::cout << "CURRENT BLOCK IS " << current_block << "\n";
        //std::cout << "UPDATING WITH " << (((current_block >> (64 - 2)) & uint64_t(3))) << "\n";
        update_rolling_hash_in(((current_block >> (64 - 2)) & uint64_t(3)));
        current_block = current_block << 2;
        bits_already_read_from_current_block += 2;
        hashed_count += 1;
        if (bits_already_read_from_current_block == 64)
        {
            current_block_position += 1;
            current_block = kmer_factory->blocks[current_block_position];
            bits_already_read_from_current_block = 0;
        }
    }
}

uint64_t RollingHasher1::update_rolling_hash(uint64_t in, uint64_t out)
{
    // Update number of characters in the hash value
    hashed_count = std::min(m, hashed_count+1);
    // Based on the outgoing character and hash count, update hash value accordingly
    if ((out == 0) || (hashed_count < m))
    {
        return update_rolling_hash_in(in);
    }
    else
    {
        return update_rolling_hash_in_and_out(in, out);
    }
}

uint64_t RollingHasher1::get_current_hash()
{
    //std::cout << "Current hash is: " << current_hash << "\n";
    return current_hash;
}

void RollingHasher1::reset()
{
    current_hash = 0;
    hashed_count = 0;
}


// Double hashing probing
uint64_t ProbeHasher1::probe_1(uint64_t item, uint64_t P)
{
    //return 1;
    return P - (item % P);
    //if (item > 0)
    //    return item % M;
    //return 1;
}

// Quadratic probing
uint64_t ProbeHasher1::probe_2(uint64_t iteration)
{
    return iteration * iteration;
    //return 1;
    //return P - (item % P);
    //if (item > 0)
    //    return item % M;
    //return 1;
}