#include "kmer_factory.hpp"

/*

    k-mer factory class implementation

*/

/*
    === Forward version starts ==========================================================================================
*/
KMerFactory2BC::KMerFactory2BC(uint64_t k)
{

    // Example of 3 blocks and how the characters fit in them
    // Leftmost block is the one with empty slots
    // (EX = Empty 2 bit slot number X)
    // (CX = Used 2 bit slot number X)

    // E01 : E02 : C01 : C02 | C03 : C04 : C05 : C06 | C07 : C08 : C09 : C10 | C11 : C12 : C13 : C14
    // C15 : C16 : C17 : C18 | C19 : C20 : C21 : C22 | C23 : C24 : C25 : C26 | C27 : C28 : C29 : C30 

    // Bits per one character
    character_bits = 2;
    // Set k-mer length
    kmer_length = k;
    // Set number of characters stored
    characters_stored = 0;
    // Set the pushed off character
    pushed_off_character = 0;
    // Set the mask fo the rightmost character of a block
    right_char_mask = uint64_t(3);
    // Set the mask fo the leftmost character of a block
    left_char_mask = right_char_mask << (64 - character_bits);
    // Calculate number of useb bits in last block
    bits_in_last_block = (character_bits*k) % 64;
    // Calculate needed blocks
    number_of_blocks = std::ceil((character_bits*k)/64.0);
    // Mask for the newest character in the rightmost block 
    right_block_right_char_mask = right_char_mask;
    // Mask for the oldest character in the leftmost block
    left_block_left_char_mask = (right_char_mask << (bits_in_last_block-character_bits));
    // Calculate mask for used bits in the leftmost block
    used_left_block_mask = 0;
    for (int b = 0; b < bits_in_last_block; b++){used_left_block_mask <<= 1; used_left_block_mask |= uint64_t(1);}
    // Initialize vector for 64 bit blocks
    blocks = new uint64_t[number_of_blocks];
    for (int i = 0; i < number_of_blocks; i++)
        blocks[i] = 0;           
}

KMerFactory2BC::~KMerFactory2BC()
{
    delete[] blocks;
}

bool KMerFactory2BC::current_kmer_is_real()
{
    return kmer_length == characters_stored;
}

int KMerFactory2BC::get_number_of_stored_characters()
{
    return characters_stored;
}

void KMerFactory2BC::reset()
{
    characters_stored = 0;
    pushed_off_character = 0;
    for (int i = 0; i < number_of_blocks; i++)
        blocks[i] = 0;
}

void KMerFactory2BC::push_new_character(char c)
{
    //std::cout << "Pushing a new character to k-mer factory: " << c << "\n" ;
    pushed_off_character = get_leftmost_character();
    //std::cout << "Pushed off character is " << pushed_off_character << "\n";
    uint64_t new_char = twobitstringfunctions::char2int(c);
    //std::cout << "To be pushed character is " << new_char << "\n";
    
    // If we are trying to push an invalid character, reset the factory
    if (new_char > uint64_t(3))
    {
        std::cout << "* Warning * Invalid character tried to get into the k-mer factory\n";
        reset();
        return;
    }

    for (int i = 0; i < number_of_blocks-1; i++)
    {
        blocks[i] <<= character_bits;
        blocks[i] |= (blocks[i+1] >> (64-character_bits)); 
    }
    blocks[number_of_blocks-1] <<= character_bits;
    blocks[number_of_blocks-1] |= new_char;
    blocks[0] &= used_left_block_mask;

    //std::cout << "Left block mask is " << used_left_block_mask << "\n";
    //std::cout << "Left block is now: " << blocks[0] << "\n";

    characters_stored = std::min(characters_stored+1, kmer_length);
}

void KMerFactory2BC::push_new_integer(uint64_t c)
{
    //std::cout << "Pushing a new character to k-mer factory: " << c << "\n" ;
    pushed_off_character = get_leftmost_character();
    //std::cout << "Pushed off character is " << pushed_off_character << "\n";
    uint64_t new_char = c;
    //std::cout << "To be pushed character is " << new_char << "\n";
    
    // If we are trying to push an invalid character, reset the factory
    if (new_char > uint64_t(3))
    {
        std::cout << "* Warning * Invalid character tried to get into the k-mer factory\n";
        reset();
        return;
    }

    for (int i = 0; i < number_of_blocks-1; i++)
    {
        blocks[i] <<= character_bits;
        blocks[i] |= (blocks[i+1] >> (64-character_bits)); 
    }
    blocks[number_of_blocks-1] <<= character_bits;
    blocks[number_of_blocks-1] |= new_char;
    blocks[0] &= used_left_block_mask;

    //std::cout << "Left block mask is " << used_left_block_mask << "\n";
    //std::cout << "Left block is now: " << blocks[0] << "\n";

    characters_stored = std::min(characters_stored+1, kmer_length);
}


uint64_t KMerFactory2BC::get_leftmost_character()
{   
    return ((blocks[0]&left_block_left_char_mask) >> (bits_in_last_block-character_bits));
}

uint64_t KMerFactory2BC::get_rightmost_character()
{
    return (blocks[number_of_blocks-1]&right_block_right_char_mask);
}

uint64_t KMerFactory2BC::get_pushed_off_character()
{
    return pushed_off_character;
}

uint64_t KMerFactory2BC::get_newest_character()
{
    return get_rightmost_character();
}

/*
    === Forward version ends ==========================================================================================
*/


/*
    === Canonical version starts ==========================================================================================
*/

KMerFactoryCanonical2BC::KMerFactoryCanonical2BC(uint64_t k)
{

    // Example of 3 blocks and how the characters fit in them
    // Leftmost block is the one with empty slots
    // (EX = Empty 2 bit slot number X)
    // (CX = Used 2 bit slot number X)

    // E01 : E02 : C01 : C02 | C03 : C04 : C05 : C06 | C07 : C08 : C09 : C10 | C11 : C12 : C13 : C14
    // C15 : C16 : C17 : C18 | C19 : C20 : C21 : C22 | C23 : C24 : C25 : C26 | C27 : C28 : C29 : C30 

    // Bits per one character
    character_bits = 2;
    // Set k-mer length
    kmer_length = k;
    // Set number of characters stored
    characters_stored = 0;
    // Set the pushed off character
    pushed_off_character_forward = 0;
    // Set the mask fo the rightmost character of a block
    right_char_mask = uint64_t(3);
    // Set the mask fo the leftmost character of a block
    left_char_mask = right_char_mask << (64 - character_bits);
    // Calculate number of useb bits in last block
    bits_in_last_block = (character_bits*k) % 64;
    // Calculate needed blocks
    number_of_blocks = std::ceil((character_bits*k)/64.0);
    // Mask for the newest character in the rightmost block 
    right_block_right_char_mask = right_char_mask;
    // Mask for the oldest character in the leftmost block
    left_block_left_char_mask = (right_char_mask << (bits_in_last_block-character_bits));
    // Calculate mask for used bits in the leftmost block
    used_left_block_mask = 0;
    for (int b = 0; b < bits_in_last_block; b++){used_left_block_mask <<= 1; used_left_block_mask |= uint64_t(1);}
    // Initialize vector for 64 bit blocks
    blocks_forward = new uint64_t[number_of_blocks];
    blocks_backward = new uint64_t[number_of_blocks];
    for (int i = 0; i < number_of_blocks; i++)
    {
        blocks_forward[i] = 0;
        blocks_backward[i] = 0;
    }
    // Set canonicality to forward
    forward_is_canonical = true;
    previous_forward_was_canonical = true;
    previous_kmer_exists = false;
                   
}

KMerFactoryCanonical2BC::~KMerFactoryCanonical2BC()
{
    delete[] blocks_forward;
    delete[] blocks_backward;
}

bool KMerFactoryCanonical2BC::current_kmer_is_real()
{
    return kmer_length == characters_stored;
}

bool KMerFactoryCanonical2BC::forward_kmer_is_canonical()
{
    return forward_is_canonical;
} 

bool KMerFactoryCanonical2BC::previous_forward_kmer_was_canonical()
{
    return previous_forward_is_canonical;
} 

bool KMerFactoryCanonical2BC::previous_kmer_existed()
{
    return previous_kmer_exists;
} 

int KMerFactoryCanonical2BC::get_number_of_stored_characters()
{
    return characters_stored;
}


void KMerFactoryCanonical2BC::reset()
{
    forward_is_canonical = true;
    previous_forward_was_canonical = true;
    previous_kmer_exists = false;
    characters_stored = 0;
    pushed_off_character = 0;
    for (int i = 0; i < number_of_blocks; i++)
    {
        blocks_forward[i] = 0;
        blocks_backward[i] = 0;
    }
}

void KMerFactoryCanonical2BC::push_new_character(char c)
{
    pushed_off_character = get_forward_leftmost_character();
    uint64_t new_char = twobitstringfunctions::char2int(c);
    uint64_t new_char_reversed = twobitstringfunctions::reverse_int(new_char);
    
    // If we are trying to push an invalid character, reset the factory
    if (new_char > uint64_t(3))
    {
        std::cout << "* Warning * Invalid character tried to get into the k-mer factory\n";
        reset();
        return;
    }

    // Forward k-mer stuff
    for (int i = 0; i < number_of_blocks-1; i++)
    {
        blocks_forward[i] <<= character_bits;
        blocks_forward[i] |= (blocks_forward[i+1] >> (64-character_bits)); 
    }
    blocks_forward[number_of_blocks-1] <<= character_bits;
    blocks_forward[number_of_blocks-1] |= new_char;
    blocks_forward[0] &=  used_left_block_mask;
    
    // Reverse k-mer stuff
    for (int i = number_of_blocks-1; i > 0; i--)
    {
        blocks_backward[i] >>= character_bits;
        blocks_backward[i] |= (blocks_backward[i-1] << (64-character_bits));
    }
    blocks_backward[0] >>= character_bits;
    blocks_backward[0] |= (new_char_reversed << (bits_in_last_block-character_bits));

    // Resolve canonical orientation
    previous_forward_was_canonical = forward_is_canonical;
    forward_is_canonical = true;
    for (int i = 0; i < number_of_blocks; i++)
    {
        if (blocks_forward[i] < blocks_backward[i])
        {
            break;
        }
        else if (blocks_forward[i] > blocks_backward[i])
        {
            forward_is_canonical = false;
            break;
        }
    }
    if (characters_stored == kmer_length)
    {
        previous_kmer_exists = true;
    }
    characters_stored = std::min(characters_stored+1, kmer_length);
}

void KMerFactoryCanonical2BC::push_new_integer(uint64_t c)
{
    pushed_off_character = get_forward_leftmost_character();
    uint64_t new_char = c;
    uint64_t new_char_reversed = twobitstringfunctions::reverse_int(new_char);
    
    // If we are trying to push an invalid character, reset the factory
    if (new_char > uint64_t(3))
    {
        std::cout << "* Warning * Invalid character tried to get into the k-mer factory\n";
        reset();
        return;
    }

    // Forward k-mer stuff
    for (int i = 0; i < number_of_blocks-1; i++)
    {
        blocks_forward[i] <<= character_bits;
        blocks_forward[i] |= (blocks_forward[i+1] >> (64-character_bits)); 
    }
    blocks_forward[number_of_blocks-1] <<= character_bits;
    blocks_forward[number_of_blocks-1] |= new_char;
    blocks_forward[0] &=  used_left_block_mask;
    
    // Reverse k-mer stuff
    for (int i = number_of_blocks-1; i > 0; i--)
    {
        blocks_backward[i] >>= character_bits;
        blocks_backward[i] |= (blocks_backward[i-1] << (64-character_bits));
    }
    blocks_backward[0] >>= character_bits;
    blocks_backward[0] |= (new_char_reversed << (bits_in_last_block-character_bits));

    // Resolve canonical orientation
    previous_forward_was_canonical = forward_is_canonical;
    forward_is_canonical = true;
    for (int i = 0; i < number_of_blocks; i++)
    {
        if (blocks_forward[i] < blocks_backward[i])
        {
            break;
        }
        else if (blocks_forward[i] > blocks_backward[i])
        {
            forward_is_canonical = false;
            break;
        }
    }
    if (characters_stored == kmer_length)
    {
        previous_kmer_exists = true;
    }
    characters_stored = std::min(characters_stored+1, kmer_length);
}


uint64_t KMerFactoryCanonical2BC::get_forward_leftmost_character()
{   
    return ((blocks_forward[0]&left_block_left_char_mask) >> (bits_in_last_block-character_bits));
}

uint64_t KMerFactoryCanonical2BC::get_forward_rightmost_character()
{
    return (blocks_forward[number_of_blocks-1]&right_block_right_char_mask);
}

uint64_t KMerFactoryCanonical2BC::get_forward_pushed_off_character()
{
    return pushed_off_character;
}

uint64_t KMerFactoryCanonical2BC::get_forward_newest_character()
{
    return get_rightmost_character();
}

uint64_t KMerFactoryCanonical2BC::get_canonical_block(uint64_t i)
{
    if (forward_is_canonical)
    {
        return blocks_forward[i];
    }
    else
    {
        return blocks_backward[i];
    }
}

uint64_t KMerFactoryCanonical2BC::get_noncanonical_block(uint64_t i)
{
    if (forward_is_canonical)
    {
        return blocks_backward[i];
    }
    else
    {
        return blocks_forward[i];
    }
}

uint64_t KMerFactoryCanonical2BC::get_canonical_right_character()
{
    return (uint64_t(3)&get_canonical_block(number_of_blocks-1));
}

uint64_t KMerFactoryCanonical2BC::get_canonical_left_character()
{
    return (uint64_t(3)&(get_canonical_block(0)>>(64-bits_in_last_block-2)));
}


/*
    === Canonical version ends ==========================================================================================
*/