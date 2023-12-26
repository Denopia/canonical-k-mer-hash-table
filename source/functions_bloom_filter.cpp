#include "functions_bloom_filter.hpp"
#include "../external/xxHash/xxhash.h"

namespace bfhf
{

    uint64_t hf_flex(uint64_t value, uint64_t size, uint64_t seed)
    {
        XXH64_hash_t h = XXH64(&value, 8, seed);
        //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf_flex_2P(uint64_t value, uint64_t size, uint64_t seed)
    {
        XXH64_hash_t h = XXH64(&value, 8, seed);
        //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h & (size-1);
    }


    uint64_t hf11(uint64_t value, uint64_t size)
    {
        uint64_t seed = 31;
        XXH64_hash_t h = XXH64(&value, 8, seed);
        //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf12(uint64_t value, uint64_t size)
    {
        uint64_t seed = 47;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf13(uint64_t value, uint64_t size)
    {
        uint64_t seed = 59;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf14(uint64_t value, uint64_t size)
    {
        uint64_t seed = 67;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf15(uint64_t value, uint64_t size)
    {
        uint64_t seed = 161;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf21(uint64_t value, uint64_t size)
    {
        uint64_t seed = 79;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }
    
    uint64_t hf22(uint64_t value, uint64_t size)
    {
        uint64_t seed = 101;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf23(uint64_t value, uint64_t size)
    {
        uint64_t seed = 127;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf24(uint64_t value, uint64_t size)
    {
        uint64_t seed = 163;
        XXH64_hash_t h = XXH64(&value, 8, seed);
	    //xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }
}