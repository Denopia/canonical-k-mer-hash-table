#include <cstdint>
#include <iostream>

#pragma once

namespace bfhf
{
    uint64_t hf_flex(uint64_t value, uint64_t size, uint64_t seed);

    uint64_t hf_flex_2P(uint64_t value, uint64_t size, uint64_t seed);
}