#include "functions_bloom_filter.hpp"

namespace bfhf
{

    uint64_t hf_flex(uint64_t value, uint64_t size, uint64_t seed)
    {
        xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf_flex_2P(uint64_t value, uint64_t size, uint64_t seed)
    {
        xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h & (size-1);
    }


    uint64_t hf11(uint64_t value, uint64_t size)
    {
        uint64_t seed = 31;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf12(uint64_t value, uint64_t size)
    {
        uint64_t seed = 47;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf13(uint64_t value, uint64_t size)
    {
        uint64_t seed = 59;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf14(uint64_t value, uint64_t size)
    {
        uint64_t seed = 67;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf15(uint64_t value, uint64_t size)
    {
        uint64_t seed = 161;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf21(uint64_t value, uint64_t size)
    {
        uint64_t seed = 79;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }
    
    uint64_t hf22(uint64_t value, uint64_t size)
    {
        uint64_t seed = 101;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf23(uint64_t value, uint64_t size)
    {
        uint64_t seed = 127;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }

    uint64_t hf24(uint64_t value, uint64_t size)
    {
        uint64_t seed = 163;
	    xxh::hash64_t h = xxh::xxhash<64>(&value, 8, seed);
        return h % size;
    }
}