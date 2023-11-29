#pragma once

#include <iostream>
#include <cstdint>
#include <cstring> // for memset



class Atomic8
{
    public:
        std::atomic<uint8_t> value;
        Atomic8():value(0){};
        ~Atomic8(){};
        void set(uint8_t nv){
            uint8_t ev = value.load(std::memory_order_acquire);
            while(!value.compare_exchange_strong(ev, nv, std::memory_order_acq_rel,std::memory_order_relaxed)){}
        }
        uint8_t get(){
            return value.load(std::memory_order_acquire);
        }
};


// Atomic bit array. INTERLEAVED AND SPLIT IN TWO
// Atomic bit array
class MyAtomicBitArrayFT
{
    private:
        std::atomic<uint8_t> * array1;
        std::atomic<uint8_t> * array2;
        uint8_t * bit_tester;
        uint64_t bits;
        uint64_t bytes;
        uint64_t half_bits;
        uint64_t half_bytes;
        bool squeezed;

    public:
        MyAtomicBitArrayFT(uint64_t size);
        ~MyAtomicBitArrayFT();
        bool test(uint64_t i);
        bool set(uint64_t i);
        void squeeze();

};

MyAtomicBitArrayFT::MyAtomicBitArrayFT(uint64_t size)
{
    squeezed = false;
    bit_tester = new uint8_t[8]{128,64,32,16,8,4,2,1};
    uint64_t divres = size >> 3;
    uint64_t rem = size & 7;
    bits = size;
    bytes = divres;
    if (rem != 0)
        bytes += 1;
    half_bits = bits >> 1;
    half_bytes = bytes >> 1;
    //bytes = (size & 7) == 0 ? (size >> 3) : (size >> 3) + 1;
    array1 = new std::atomic<uint8_t>[bytes>>1]{};
    array2 = new std::atomic<uint8_t>[bytes>>1]{};
}

MyAtomicBitArrayFT::~MyAtomicBitArrayFT()
{
    if (!squeezed)
        delete[] array2;
    delete[] array1;
    delete[] bit_tester;
}

// Test if bit is set
bool MyAtomicBitArrayFT::test(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;
    if (divres < half_bytes)
        return (bit_tester[rem]) == (bit_tester[rem] & array1[divres].load(std::memory_order_acquire));
    else
        return (bit_tester[rem]) == (bit_tester[rem] & array2[divres-half_bytes].load(std::memory_order_acquire));
}

// Set bit to 1. Return true if bit was set, return false if bit could not be set (was set by someone else).
bool MyAtomicBitArrayFT::set(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;

    if (divres < half_bytes){
        uint8_t expected_data = array1[divres].load(std::memory_order_acquire);
        while(!array1[divres].compare_exchange_strong(
            expected_data, 
            expected_data | bit_tester[rem],
            std::memory_order_acq_rel,
            std::memory_order_relaxed))
        {
            // If bit is already set, return "failed to set"
            if ((expected_data&bit_tester[rem]) == bit_tester[rem])
            {
                //std::cout << "SOMEONE ALREADY FLIPPED IT\n";
                return false;
            }
        }
    } else {
        uint8_t expected_data = array2[divres-half_bytes].load(std::memory_order_acquire);
        while(!array2[divres-half_bytes].compare_exchange_strong(
            expected_data, 
            expected_data | bit_tester[rem],
            std::memory_order_acq_rel,
            std::memory_order_relaxed))
        {
            // If bit is already set, return "failed to set"
            if ((expected_data&bit_tester[rem]) == bit_tester[rem])
            {
                //std::cout << "SOMEONE ALREADY FLIPPED IT\n";
                return false;
            }
        }
    }
    // If we get here, return "bit set succesful"
    return true;
}


void MyAtomicBitArrayFT::squeeze()
{
    squeezed = true;
    for (int i = 0; 2*i+1 < bytes; i++){
        uint8_t new_value = 0;

        if (2*i < half_bytes){
            new_value = new_value | ((array1[2*i].load(std::memory_order_acquire)&64)<<1);
            new_value = new_value | ((array1[2*i].load(std::memory_order_acquire)&16)<<2);
            new_value = new_value | ((array1[2*i].load(std::memory_order_acquire)&4)<<3);
            new_value = new_value | ((array1[2*i].load(std::memory_order_acquire)&1)<<4);
        } else {
            new_value = new_value | ((array2[2*i - half_bytes].load(std::memory_order_acquire)&64)<<1);
            new_value = new_value | ((array2[2*i - half_bytes].load(std::memory_order_acquire)&16)<<2);
            new_value = new_value | ((array2[2*i - half_bytes].load(std::memory_order_acquire)&4)<<3);
            new_value = new_value | ((array2[2*i - half_bytes].load(std::memory_order_acquire)&1)<<4);
        }
        
         if (2*i < half_bytes){
            new_value = new_value | ((array1[2*i+1].load(std::memory_order_acquire)&64)>>3);
            new_value = new_value | ((array1[2*i+1].load(std::memory_order_acquire)&16)>>2);
            new_value = new_value | ((array1[2*i+1].load(std::memory_order_acquire)&4)>>1);
            new_value = new_value | ((array1[2*i+1].load(std::memory_order_acquire)&1));
        } else {
            new_value = new_value | ((array2[2*i+1 - half_bytes].load(std::memory_order_acquire)&64)>>3);
            new_value = new_value | ((array2[2*i+1 - half_bytes].load(std::memory_order_acquire)&16)>>2);
            new_value = new_value | ((array2[2*i+1 - half_bytes].load(std::memory_order_acquire)&4)>>1);
            new_value = new_value | ((array2[2*i+1 - half_bytes].load(std::memory_order_acquire)&1));
        } 

        uint8_t old_value = array1[i].load(std::memory_order_acquire);

        while(!array1[i].compare_exchange_strong(old_value, new_value, std::memory_order_acq_rel, std::memory_order_relaxed)){}
    }
    delete[] array2;
}


// Atomic bit array
class MyAtomicBitVector
{
    private:
        std::vector<Atomic8> bitvector;
        uint8_t * bit_tester;
        uint64_t bits;
        uint64_t bytes;
        bool squeezed;
        uint64_t size;

    public:
        MyAtomicBitVector(uint64_t s);
        ~MyAtomicBitVector();
        bool test(uint64_t i);
        bool set(uint64_t i);
        void squeeze();

};

MyAtomicBitVector::MyAtomicBitVector(uint64_t s)
{
    size = s;
    squeezed = false;
    bit_tester = new uint8_t[8]{128,64,32,16,8,4,2,1};
    uint64_t divres = size >> 3;
    uint64_t rem = size & 7;
    bits = size;
    bytes = divres;
    if (rem != 0)
        bytes += 1;
    bitvector = std::vector<Atomic8>(bytes);
}

MyAtomicBitVector::~MyAtomicBitVector()
{
    delete[] bit_tester;
}

// Test if bit is set
bool MyAtomicBitVector::test(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;
    return (bit_tester[rem]) == (bit_tester[rem] & bitvector[divres].value.load(std::memory_order_acquire));
}

// Set bit to 1. Return true if bit was set, return false if bit could not be set (was set by someone else).
bool MyAtomicBitVector::set(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;
    uint8_t expected_data = bitvector[divres].value.load(std::memory_order_acquire);
    
    while(!bitvector[divres].value.compare_exchange_strong(
        expected_data, 
        expected_data | bit_tester[rem],
        std::memory_order_acq_rel,
        std::memory_order_relaxed))
    {
        // If bit is already set, return "failed to set"
        if ((expected_data&bit_tester[rem]) == bit_tester[rem])
        {
            //std::cout << "SOMEONE ALREADY FLIPPED IT\n";
            return false;
        }
    }
    // If we get here, return "bit set succesful"
    return true;
}

void MyAtomicBitVector::squeeze()
{
    squeezed = true;
    for (int i = 0; 2*i+1 < bytes; i++){
        uint8_t new_value = 0;

        new_value = new_value | ((bitvector[2*i].value.load(std::memory_order_acquire)&64)<<1);
        new_value = new_value | ((bitvector[2*i].value.load(std::memory_order_acquire)&16)<<2);
        new_value = new_value | ((bitvector[2*i].value.load(std::memory_order_acquire)&4)<<3);
        new_value = new_value | ((bitvector[2*i].value.load(std::memory_order_acquire)&1)<<4);

        new_value = new_value | ((bitvector[2*i+1].value.load(std::memory_order_acquire)&64)>>3);
        new_value = new_value | ((bitvector[2*i+1].value.load(std::memory_order_acquire)&16)>>2);
        new_value = new_value | ((bitvector[2*i+1].value.load(std::memory_order_acquire)&4)>>1);
        new_value = new_value | ((bitvector[2*i+1].value.load(std::memory_order_acquire)&1));

        uint8_t old_value = bitvector[i].value.load(std::memory_order_acquire);

        while(!bitvector[i].value.compare_exchange_strong(old_value, new_value, std::memory_order_acq_rel, std::memory_order_relaxed)){}
    }
    bytes = bytes / 2;
    for (int i = 0; i < bytes; i++)
        bitvector.pop_back();
    bitvector.shrink_to_fit();
}


// Atomic bit array
class MyAtomicBitArray
{
    private:
        std::atomic<uint8_t> * array;
        uint8_t * bit_tester;
        uint64_t bits;
        uint64_t bytes;

    public:
        MyAtomicBitArray(uint64_t size);
        ~MyAtomicBitArray();
        bool test(uint64_t i);
        bool set(uint64_t i);

};

MyAtomicBitArray::MyAtomicBitArray(uint64_t size)
{
    bit_tester = new uint8_t[8]{128,64,32,16,8,4,2,1};
    uint64_t divres = size >> 3;
    uint64_t rem = size & 7;
    bits = size;
    bytes = divres;
    if (rem != 0)
        bytes += 1;
    //bytes = (size & 7) == 0 ? (size >> 3) : (size >> 3) + 1;
    array = new std::atomic<uint8_t>[bytes]{};
}

MyAtomicBitArray::~MyAtomicBitArray()
{
    delete[] array;
    delete[] bit_tester;
}

// Test if bit is set
bool MyAtomicBitArray::test(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;
    return (bit_tester[rem]) == (bit_tester[rem] & array[divres].load(std::memory_order_acquire));
}

// Set bit to 1. Return true if bit was set, return false if bit could not be set (was set by someone else).
bool MyAtomicBitArray::set(uint64_t i)
{
    uint64_t divres = i >> 3;
    uint64_t rem = i & 7;
    uint8_t expected_data = array[divres].load(std::memory_order_acquire);
    
    while(!array[divres].compare_exchange_strong(
        expected_data, 
        expected_data | bit_tester[rem],
        std::memory_order_acq_rel,
        std::memory_order_relaxed))
    {
        // If bit is already set, return "failed to set"
        if ((expected_data&bit_tester[rem]) == bit_tester[rem])
        {
            //std::cout << "SOMEONE ALREADY FLIPPED IT\n";
            return false;
        }
    }
    // If we get here, return "bit set succesful"
    return true;
}




// Regular bit array
class MyBitArray
{
    private:
        uint16_t * array;
        uint16_t * bit_tester;
        uint64_t bits;
        uint64_t bytes;

    public:
        MyBitArray(uint64_t size);
        ~MyBitArray();
        bool test(uint64_t i);
        void set(uint64_t i);

};

MyBitArray::MyBitArray(uint64_t size)
{
    //bit_tester = new uint8_t[8]{128,64,32,16,8,4,2,1};
    bit_tester = new uint16_t[16]{1<<15,1<<14,1<<13,1<<12,1<<11,1<<10,1<<9,1<<8,1<<7,1<<6,1<<5,1<<4,1<<3,1<<2,1<<1,1};
    uint64_t divres = size >> 4; // 3
    uint64_t rem = size & 15; // 7
    bits = size;
    bytes = divres;
    if (rem != 0)
        bytes += 1;
    //bytes = (size & 7) == 0 ? (size >> 3) : (size >> 3) + 1;
    array = new uint16_t[bytes]{};
}

MyBitArray::~MyBitArray()
{
    delete[] array;
    delete[] bit_tester;
}

bool MyBitArray::test(uint64_t i)
{
    uint64_t divres = i >> 4;
    uint64_t rem = i & 15;
    return (bit_tester[rem]) == (bit_tester[rem]&array[divres]);
}

void MyBitArray::set(uint64_t i)
{
    uint64_t divres = i >> 4;
    uint64_t rem = i & 15;
    array[divres] |= bit_tester[rem]; 
}
