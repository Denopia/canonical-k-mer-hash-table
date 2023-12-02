#include <iostream>
#include <cstdint>
#include <xxhash.h>
#include <vector>
//#include <sdsl/bit_vectors.hpp>
#include <random>  // Include the <random> header
#include "mybitarray.hpp"


#pragma once

// === SIMPLE RESIZING ===
// Atomic double bloom filter using my atomic bitarray
class DoubleAtomicDoubleBloomFilterS {
private:
    bool resized;
    std::size_t mask;
    std::size_t numHashFunctions;
    std::atomic<std::size_t> new_in_first;
    std::atomic<std::size_t> new_in_second;
    std::vector<std::size_t> seeds;  // Store random seeds
    std::atomic<std::size_t> failed_insertions_in_first;
    MyAtomicBitArray * bitArray2;
    MyAtomicBitArray * bitArray1;
    

public:
    
    DoubleAtomicDoubleBloomFilterS(std::size_t size, std::size_t numHashFunctions)
        : resized(false), mask(size - 1), numHashFunctions(numHashFunctions), new_in_first(0), new_in_second(0), failed_insertions_in_first(0) {
        seeds = generate_seeds_2(numHashFunctions);  // Generate random seeds
        bitArray2 = new MyAtomicBitArray(size);
        bitArray1 = new MyAtomicBitArray(size);
        //std::cout << "numver of hash functionsis " << numHashFunctions << "\n";
    }

    ~DoubleAtomicDoubleBloomFilterS(){
        delete bitArray2;
        if (!resized)
            delete bitArray1;
    }

    uint64_t get_failed_insertions_in_first(){
        return failed_insertions_in_first.load(std::memory_order_acquire);
    }

    uint64_t get_new_in_first(){
        return new_in_first.load(std::memory_order_acquire);
    }

    uint64_t get_new_in_second(){
        return new_in_second.load(std::memory_order_acquire);
    }

    inline void calculate_hashes(uint64_t value, std::vector<uint64_t> &hash_values){
        for (std::size_t i = 0; i < numHashFunctions; ++i) {
            std::size_t hash = XXH64(&value, sizeof(value), seeds[i]);
            hash_values[i] = hash & mask;
        }
    }

    void resize(){
        resized = true;
        delete bitArray1;
    }

    uint64_t first_contains(std::vector<uint64_t> &hash_values){
        if (!resized){
            uint64_t set_bits = 0;
            for (auto hv : hash_values) {
                if(bitArray1->test(hv)){
                    set_bits+=1;
                }
            }
            return set_bits;
        } else {
            // First bloom filter does not exist anymore
            return false;
        }
        
    }

    uint64_t second_contains(std::vector<uint64_t> &hash_values){
        uint64_t set_bits = 0;
        for (auto hv : hash_values) {
            if(bitArray2->test(hv)){
                set_bits+=1;
            }
        }
        return set_bits;
    }

    bool insert_in_first(std::vector<uint64_t> &hash_values, uint64_t already_set_bits){
        //std::cout << "== Inserting in first\n";
        uint64_t my_set_bits = 0;
        //std::cout << "== Number of bits I have set: " << my_set_bits << "\n";
        for (auto hv : hash_values) {
            //std::cout << "== Inserting bit " << hv << "\n";
            if(!bitArray1->test(hv)){
                //std::cout << "--- It does not exists, so I am inserting it\n";
                if(bitArray1->set(hv)){
                    my_set_bits += 1;
                    //std::cout << "--- I inserted it succesfully and Increase my set bit count to " << my_set_bits << "\n";
                } else {
                    //std::cout << "--- Inserting failed, I cannot Increase my set bits count\n";
                }
            }
        }
        //std::cout << "== All bits have been set.\n";
        //std::cout << "== Number of hash functions is " << numHashFunctions << "\n";
        //std::cout << "== Already set bits was " << already_set_bits << "\n";
        //std::cout << "== My set bits is " << my_set_bits << "\n";
        //if (my_set_bits == (numHashFunctions-already_set_bits))
        //    std::cout << "== Insertion result is: success\n";
        //else
        //    std::cout << "== Insertion result is: fail\n"; 
        return (my_set_bits == (numHashFunctions-already_set_bits));
    }

    bool insert_in_second(std::vector<uint64_t> &hash_values, uint64_t already_set_bits){
        uint64_t my_set_bits = 0;
        for (auto hv : hash_values) {
            if(!bitArray2->test(hv)){
                if(bitArray2->set(hv)){
                    my_set_bits += 1;
                }
            }
        }
        return (my_set_bits == (numHashFunctions-already_set_bits));
    }

    // Needs to be modified to take atomicity into account
    void insertion_process(uint64_t value, std::vector<uint64_t> &hash_values){
        
        // Calculate hashes
        calculate_hashes(value, hash_values); 
        //std::cout << "HASHES ARE\n";
        //for (auto a : hash_values)
        //    std::cout << a << "\n";
        // Check if exists in second first
        uint64_t set_second_bits = second_contains(hash_values);
        if (set_second_bits == numHashFunctions){
            //std::cout << "Found in second1\n";
            return;
        }
        // Check if exists in first
        uint64_t set_first_bits = first_contains(hash_values);
        //std::cout << "!!BEFORE!! CONTAINED " << set_first_bits << " BITS IN FIRST AND " << set_second_bits << " IN SECOND\n";
        if (set_first_bits == numHashFunctions){
            //std::cout << "Found in first1\n";
            if (insert_in_second(hash_values, set_second_bits)){
                // If insertion was succesful, increase count
                //std::cout << "Increase second count1\n";
                uint64_t nisu = new_in_second.load(std::memory_order_acquire);
                while(!new_in_second.compare_exchange_strong(nisu, nisu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
            }
        } else {
            // If not, insert in first
            //std::cout << "inserting in first\n";
            if (insert_in_first(hash_values, set_first_bits)){
                //std::cout << "-inserted succesfully\n";
                // If insertion was succesful, increase count
                uint64_t nifu = new_in_first.load(std::memory_order_acquire);
                while(!new_in_first.compare_exchange_strong(nifu, nifu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
            // It was inserted by someone else so we need to put it in the second filter
            } else {
                uint64_t fails = failed_insertions_in_first.load(std::memory_order_acquire);
                while(!failed_insertions_in_first.compare_exchange_strong(fails, fails+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
                //std::cout << "-insertion failed\n";
                if (insert_in_second(hash_values, set_second_bits)){
                    // If insertion was succesful, increase count
                    //std::cout << "Increase second count2\n";
                    uint64_t nisu = new_in_second.load(std::memory_order_acquire);
                    while(!new_in_second.compare_exchange_strong(nisu, nisu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
                }
            }
        }
        //set_first_bits = first_contains(hash_values);
        //set_second_bits = second_contains(hash_values);
        //std::cout << "!!AFTER!! CONTAINED " << set_first_bits << " BITS IN FIRST AND " << set_second_bits << " IN SECOND\n";
    }

    

private:
    // Function to generate random seeds
    std::vector<std::size_t> generate_seeds(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<std::size_t> dis_seed;

        for (std::size_t i = 0; i < numSeeds; ++i) {
            result.push_back(dis_seed(gen));
        }

        return result;
    }

    std::vector<std::size_t> generate_seeds_2(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::vector<uint64_t> seeds = {2411, 3253, 1061, 1129, 2269, 7309, 3491, 8237, 6359, 8779,
                                        6553, 5443, 2447, 8999, 8623, 5779, 1879, 2357, 5087, 5393,
                                        2203, 8597, 8629, 7727, 2819, 1789, 7757, 6079, 9371, 2957,
                                        2389, 4133, 4931, 2083, 8291, 1151, 4759, 7649, 6803, 1753,
                                        9613, 1979, 1877, 5479, 4799, 5303, 1759, 4451, 7841, 3461,
                                        2207, 1289, 5233, 6823, 7043, 3251, 8039, 4519, 2551, 5693,
                                        7681, 1607, 4679, 4729, 8231, 4139, 7457, 6221, 2377, 7151,
                                        3083, 6947, 7331, 3947, 6011, 7753, 2843, 3191, 7993, 4943,
                                        5801, 9901, 4001, 1933, 7523, 5273, 1721, 1093, 2579, 9719,
                                        4481, 2417, 9341, 9137, 5113, 9719, 5399, 5231, 1979, 6701,
                                        4133, 1723, 1931, 3257, 8861, 8539, 4877, 2207, 7151, 5279};
        uint64_t start = 0;
        while (start < numSeeds){
            result.push_back(seeds[start]);
            start+=1;
        }
        return result;
    }
};



// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// THIS IS THE WORKING VERSION, THE OTHER ONE IS USED FOR TESTING SOMETHING
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Atomic double bloom filter using my atomic bitarray
class DoubleAtomicDoubleBloomFilter {
private:
    //MyAtomicBitArray * bitArray;
    MyAtomicBitArrayFT * bitArray;
    std::size_t mask;
    std::size_t numHashFunctions;
    std::atomic<std::size_t> new_in_first;
    std::atomic<std::size_t> new_in_second;
    std::vector<std::size_t> seeds;  // Store random seeds
    std::atomic<std::size_t> failed_insertions_in_first;
    bool resized;

public:
    /*
    DoubleAtomicDoubleBloomFilter(std::size_t size, std::size_t numHashFunctions)
        : bitArray(2 * size), mask(size - 1), numHashFunctions(numHashFunctions), new_in_first(0), new_in_second(0), failed_insertions_in_first(0) {
        seeds = generate_seeds_2(numHashFunctions);  // Generate random seeds
        //std::cout << "numver of hash functionsis " << numHashFunctions << "\n";
    }
    */
    DoubleAtomicDoubleBloomFilter(std::size_t size, std::size_t numHashFunctions)
        : resized(false), mask(size - 1), numHashFunctions(numHashFunctions), new_in_first(0), new_in_second(0), failed_insertions_in_first(0) {
        seeds = generate_seeds_2(numHashFunctions);  // Generate random seeds
        bitArray = new MyAtomicBitArrayFT(2*size);
        //std::cout << "numver of hash functionsis " << numHashFunctions << "\n";
    }

    ~DoubleAtomicDoubleBloomFilter(){
        delete bitArray;
    }

    uint64_t get_failed_insertions_in_first(){
        return failed_insertions_in_first.load(std::memory_order_acquire);
    }

    uint64_t get_new_in_first(){
        return new_in_first.load(std::memory_order_acquire);
    }

    uint64_t get_new_in_second(){
        return new_in_second.load(std::memory_order_acquire);
    }

    inline void calculate_hashes(uint64_t value, std::vector<uint64_t> &hash_values){
        for (std::size_t i = 0; i < numHashFunctions; ++i) {
            std::size_t hash = XXH64(&value, sizeof(value), seeds[i]);
            hash_values[i] = hash & mask;
        }
    }

    /*
    void resize_OLD(){
        MyAtomicBitArray * newBitArray = new MyAtomicBitArray(mask+1);
        for (int i = 0; i < mask+1; i++){
            if (bitArray->test(2*i+1)){
                newBitArray->set(i);
            }
        }
        delete bitArray;
        bitArray = newBitArray;
        resized = true;
    }
    */
   
    void resize(){
        bitArray->squeeze();
        resized = true;
    }
    

    uint64_t first_contains(std::vector<uint64_t> &hash_values){
        if (!resized){
            uint64_t set_bits = 0;
            for (auto hv : hash_values) {
                if(bitArray->test(2*hv)){
                    set_bits+=1;
                }
            }
            return set_bits;
        } else {
            // First bloom filter does not exist anymore
            return false;
        }
        
    }

    uint64_t second_contains(std::vector<uint64_t> &hash_values){
        if (!resized){
            uint64_t set_bits = 0;
            for (auto hv : hash_values) {
                if(bitArray->test(2*hv+1)){
                    set_bits+=1;
                }
            }
            return set_bits;
        } else {
            uint64_t set_bits = 0;
            for (auto hv : hash_values) {
                if(bitArray->test(hv)){
                    set_bits+=1;
                }
            }
            return set_bits;
        }
    }

    bool insert_in_first(std::vector<uint64_t> &hash_values, uint64_t already_set_bits){
        //std::cout << "== Inserting in first\n";
        uint64_t my_set_bits = 0;
        //std::cout << "== Number of bits I have set: " << my_set_bits << "\n";
        for (auto hv : hash_values) {
            //std::cout << "== Inserting bit " << hv << "\n";
            if(!bitArray->test(2*hv)){
                //std::cout << "--- It does not exists, so I am inserting it\n";
                if(bitArray->set(2*hv)){
                    my_set_bits += 1;
                    //std::cout << "--- I inserted it succesfully and Increase my set bit count to " << my_set_bits << "\n";
                } else {
                    //std::cout << "--- Inserting failed, I cannot Increase my set bits count\n";
                }
            }
        }
        return (my_set_bits == (numHashFunctions-already_set_bits));
    }

    bool insert_in_second(std::vector<uint64_t> &hash_values, uint64_t already_set_bits){
        uint64_t my_set_bits = 0;
        for (auto hv : hash_values) {
            if(!bitArray->test(2*hv+1)){
                if(bitArray->set(2*hv+1)){
                    my_set_bits += 1;
                }
            }
        }
        return (my_set_bits == (numHashFunctions-already_set_bits));
    }

    // Needs to be modified to take atomicity into account
    void insertion_process(uint64_t value, std::vector<uint64_t> &hash_values){
        
        // Calculate hashes
        calculate_hashes(value, hash_values); 
        // Check if exists in second first
        uint64_t set_second_bits = second_contains(hash_values);
        if (set_second_bits == numHashFunctions){
            //std::cout << "Found in second1\n";
            return;
        }
        // Check if exists in first
        uint64_t set_first_bits = first_contains(hash_values);
        //std::cout << "!!BEFORE!! CONTAINED " << set_first_bits << " BITS IN FIRST AND " << set_second_bits << " IN SECOND\n";
        if (set_first_bits == numHashFunctions){
            //std::cout << "Found in first1\n";
            if (insert_in_second(hash_values, set_second_bits)){
                // If insertion was succesful, increase count
                //std::cout << "Increase second count1\n";
                uint64_t nisu = new_in_second.load(std::memory_order_acquire);
                while(!new_in_second.compare_exchange_strong(nisu, nisu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
            }
        } else {
            // If not, insert in first
            //std::cout << "inserting in first\n";
            if (insert_in_first(hash_values, set_first_bits)){
                //std::cout << "-inserted succesfully\n";
                // If insertion was succesful, increase count
                uint64_t nifu = new_in_first.load(std::memory_order_acquire);
                while(!new_in_first.compare_exchange_strong(nifu, nifu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
            // It was inserted by someone else so we need to put it in the second filter
            } else {
                uint64_t fails = failed_insertions_in_first.load(std::memory_order_acquire);
                while(!failed_insertions_in_first.compare_exchange_strong(fails, fails+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
                //std::cout << "-insertion failed\n";
                if (insert_in_second(hash_values, set_second_bits)){
                    // If insertion was succesful, increase count
                    //std::cout << "Increase second count2\n";
                    uint64_t nisu = new_in_second.load(std::memory_order_acquire);
                    while(!new_in_second.compare_exchange_strong(nisu, nisu+1, std::memory_order_acq_rel,std::memory_order_relaxed)){}
                }
            }
        }
    }

    

private:
    // Function to generate random seeds
    std::vector<std::size_t> generate_seeds(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<std::size_t> dis_seed;

        for (std::size_t i = 0; i < numSeeds; ++i) {
            result.push_back(dis_seed(gen));
        }

        return result;
    }

    std::vector<std::size_t> generate_seeds_2(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::vector<uint64_t> seeds = {2411, 3253, 1061, 1129, 2269, 7309, 3491, 8237, 6359, 8779,
                                        6553, 5443, 2447, 8999, 8623, 5779, 1879, 2357, 5087, 5393,
                                        2203, 8597, 8629, 7727, 2819, 1789, 7757, 6079, 9371, 2957,
                                        2389, 4133, 4931, 2083, 8291, 1151, 4759, 7649, 6803, 1753,
                                        9613, 1979, 1877, 5479, 4799, 5303, 1759, 4451, 7841, 3461,
                                        2207, 1289, 5233, 6823, 7043, 3251, 8039, 4519, 2551, 5693,
                                        7681, 1607, 4679, 4729, 8231, 4139, 7457, 6221, 2377, 7151,
                                        3083, 6947, 7331, 3947, 6011, 7753, 2843, 3191, 7993, 4943,
                                        5801, 9901, 4001, 1933, 7523, 5273, 1721, 1093, 2579, 9719,
                                        4481, 2417, 9341, 9137, 5113, 9719, 5399, 5231, 1979, 6701,
                                        4133, 1723, 1931, 3257, 8861, 8539, 4877, 2207, 7151, 5279};
        uint64_t start = 0;
        while (start < numSeeds){
            result.push_back(seeds[start]);
            start+=1;
        }
        return result;
    }
};

// Regular double bloom filter using my bit array
class DoubleDoubleBloomFilter {
private:
    //sdsl::bit_vector bitArray;
    MyBitArray bitArray;
    std::size_t mask;
    std::size_t numHashFunctions;
    std::size_t new_in_first;
    std::size_t new_in_second;
    std::vector<std::size_t> seeds;  // Store random seeds

public:
    DoubleDoubleBloomFilter(std::size_t size, std::size_t numHashFunctions)
        : bitArray(2 * size), mask(size - 1), numHashFunctions(numHashFunctions) {
        seeds = generate_seeds_2(numHashFunctions);  // Generate random seeds
    }

    uint64_t get_new_in_first(){
        return new_in_first;
    }

    uint64_t get_new_in_second(){
        return new_in_second;
    }

    inline void calculate_hashes(uint64_t value, std::vector<uint64_t> &hash_values){
        for (std::size_t i = 0; i < numHashFunctions; ++i) {
            std::size_t hash = XXH64(&value, sizeof(value), seeds[i]);
            hash_values[i] = hash & mask;
        }
    }

    bool first_contains(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv)){
                return false;
            }
        }
        return true;
    }

    bool second_contains(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv+1)){
                return false;
            }
        }
        return true;
    }

    bool first_contains_and_set(std::vector<uint64_t> &hash_values){
        bool found = true;
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv)){
                found = false;
                bitArray.set(2*hv);
            }
        }
        return found;
    }

    bool second_contains_and_set(std::vector<uint64_t> &hash_values){
        bool found = true;
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv+1)){
                found = false;
                bitArray.set(2*hv+1);
            }
        }
        return found;
    }

    void insert_in_first(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv)){
                bitArray.set(2*hv);
            }
        }
    }

    void insert_in_second(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray.test(2*hv+1)){
                bitArray.set(2*hv+1);
            }
        }
    }

    void insertion_process(uint64_t value, std::vector<uint64_t> &hash_values){
        
        // Calculate hashes
        calculate_hashes(value, hash_values); 
        // Check if exists in second first
        if (second_contains(hash_values)){
            return;
        }
        // Check if exists in first and also set the bits
        if (first_contains_and_set(hash_values)){
            // If already exists, insert in second
            insert_in_second(hash_values);
            new_in_second += 1;
        } else {
            new_in_first += 1;
        }
    }

private:
    // Function to generate random seeds
    std::vector<std::size_t> generate_seeds(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<std::size_t> dis_seed;

        for (std::size_t i = 0; i < numSeeds; ++i) {
            result.push_back(dis_seed(gen));
        }

        return result;
    }

    std::vector<std::size_t> generate_seeds_2(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::vector<uint64_t> seeds = {2411, 3253, 1061, 1129, 2269, 7309, 3491, 8237, 6359, 8779,
                                        6553, 5443, 2447, 8999, 8623, 5779, 1879, 2357, 5087, 5393,
                                        2203, 8597, 8629, 7727, 2819, 1789, 7757, 6079, 9371, 2957,
                                        2389, 4133, 4931, 2083, 8291, 1151, 4759, 7649, 6803, 1753,
                                        9613, 1979, 1877, 5479, 4799, 5303, 1759, 4451, 7841, 3461,
                                        2207, 1289, 5233, 6823, 7043, 3251, 8039, 4519, 2551, 5693,
                                        7681, 1607, 4679, 4729, 8231, 4139, 7457, 6221, 2377, 7151,
                                        3083, 6947, 7331, 3947, 6011, 7753, 2843, 3191, 7993, 4943,
                                        5801, 9901, 4001, 1933, 7523, 5273, 1721, 1093, 2579, 9719,
                                        4481, 2417, 9341, 9137, 5113, 9719, 5399, 5231, 1979, 6701,
                                        4133, 1723, 1931, 3257, 8861, 8539, 4877, 2207, 7151, 5279};
        uint64_t start = 0;
        while (start < numSeeds){
            result.push_back(seeds[start]);
            start+=1;
        }
        return result;
    }
};


/*
class DoubleDoubleBloomFilterS {
private:
    sdsl::bit_vector bitArray;
    std::size_t mask;
    std::size_t numHashFunctions;
    std::size_t new_in_first;
    std::size_t new_in_second;
    std::vector<std::size_t> seeds;  // Store random seeds

public:
    DoubleDoubleBloomFilterS(std::size_t size, std::size_t numHashFunctions)
        : bitArray(2 * size, 0), mask(size - 1), numHashFunctions(numHashFunctions) {
        seeds = generate_seeds_2(numHashFunctions);  // Generate random seeds
    }

    uint64_t get_new_in_first(){
        return new_in_first;
    }

    uint64_t get_new_in_second(){
        return new_in_second;
    }

    inline void calculate_hashes(uint64_t value, std::vector<uint64_t> &hash_values){
        for (std::size_t i = 0; i < numHashFunctions; ++i) {
            std::size_t hash = XXH64(&value, sizeof(value), seeds[i]);
            hash_values[i] = hash & mask;
        }
    }

    bool first_contains(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray[2*hv]){
                return false;
            }
        }
        return true;
    }

    bool second_contains(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray[2*hv+1]){
                return false;
            }
        }
        return true;
    }

    bool first_contains_and_set(std::vector<uint64_t> &hash_values){
        bool found = true;
        for (auto hv : hash_values) {
            if(!bitArray[2*hv]){
                found = false;
                bitArray[2*hv] = 1;
            }
        }
        return found;
    }

    bool second_contains_and_set(std::vector<uint64_t> &hash_values){
        bool found = true;
        for (auto hv : hash_values) {
            if(!bitArray[2*hv+1]){
                found = false;
                bitArray[2*hv+1] = 1;
            }
        }
        return found;
    }

    void insert_in_first(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray[2*hv]){
                bitArray[2*hv] = 1;
            }
        }
    }

    void insert_in_second(std::vector<uint64_t> &hash_values){
        for (auto hv : hash_values) {
            if(!bitArray[2*hv+1]){
                bitArray[2*hv+1] = 1;
            }
        }
    }

    void insertion_process(uint64_t value, std::vector<uint64_t> &hash_values){
        
        // Calculate hashes
        calculate_hashes(value, hash_values); 
        // Check if exists in second first
        if (second_contains(hash_values)){
            return;
        }
        // Check if exists in first and also set the bits
        if (first_contains_and_set(hash_values)){
            // If already exists, insert in second
            insert_in_second(hash_values);
            new_in_second += 1;
        } else {
            new_in_first += 1;
        }
    }

private:
    // Function to generate random seeds
    std::vector<std::size_t> generate_seeds(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<std::size_t> dis_seed;

        for (std::size_t i = 0; i < numSeeds; ++i) {
            result.push_back(dis_seed(gen));
        }

        return result;
    }

    std::vector<std::size_t> generate_seeds_2(std::size_t numSeeds) {
        std::vector<std::size_t> result;
        std::vector<uint64_t> seeds = {2411, 3253, 1061, 1129, 2269, 7309, 3491, 8237, 6359, 8779,
                                        6553, 5443, 2447, 8999, 8623, 5779, 1879, 2357, 5087, 5393,
                                        2203, 8597, 8629, 7727, 2819, 1789, 7757, 6079, 9371, 2957,
                                        2389, 4133, 4931, 2083, 8291, 1151, 4759, 7649, 6803, 1753,
                                        9613, 1979, 1877, 5479, 4799, 5303, 1759, 4451, 7841, 3461,
                                        2207, 1289, 5233, 6823, 7043, 3251, 8039, 4519, 2551, 5693,
                                        7681, 1607, 4679, 4729, 8231, 4139, 7457, 6221, 2377, 7151,
                                        3083, 6947, 7331, 3947, 6011, 7753, 2843, 3191, 7993, 4943,
                                        5801, 9901, 4001, 1933, 7523, 5273, 1721, 1093, 2579, 9719,
                                        4481, 2417, 9341, 9137, 5113, 9719, 5399, 5231, 1979, 6701,
                                        4133, 1723, 1931, 3257, 8861, 8539, 4877, 2207, 7151, 5279};
        int start = 0;
        while (start < numSeeds){
            result.push_back(seeds[start]);
            start+=1;
        }
        return result;
    }
};
*/

