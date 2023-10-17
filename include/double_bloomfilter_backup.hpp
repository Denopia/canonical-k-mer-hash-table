#include <iostream>
#include <cstdint>
#include <xxhash.h>
#include <vector>
//#include <sdsl/bit_vectors.hpp>
#include <random>  // Include the <random> header

#pragma once

class GPTDoubleBloomFilter {
private:
    sdsl::bit_vector bitArray;
    std::size_t mask;
    std::size_t numHashFunctions;
    std::size_t new_in_first;
    std::size_t new_in_second;
    std::vector<std::size_t> seeds;  // Store random seeds
    bool previous_in_second;

public:
    GPTDoubleBloomFilter(std::size_t size, std::size_t numHashFunctions)
        : bitArray(2 * size, 0), mask(size - 1), numHashFunctions(numHashFunctions), previous_in_second(false) {
        seeds = generate_seeds(numHashFunctions);  // Generate random seeds
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
};
