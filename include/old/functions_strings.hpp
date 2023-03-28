#include <iostream>
#include <cstdint>

#pragma once

/*

    File containing some string functions,
    currently mainly conversions between
    string and integer forms

*/

// Pure string functions
namespace purestringfunctions
{
    char reverse_char(char c);

    void reverse_this_string(std::string & s);

    int is_canonical(std::string & s);

}



// 2 bit character functions
namespace twobitstringfunctions
{

    uint64_t char2int(char c);

    char int2char(uint64_t c);

    uint64_t reverse_int(uint64_t c);

    std::string int2string_single(uint64_t s, int l);

    std::string int2string_multi(uint64_t* d, int b, int l);
    
    bool is_clean_string(std::string s);

}