#include "functions_math.hpp"


namespace mathfunctions
{
    
    uint64_t next_prime(uint64_t at_least)
    {
        // The prime candidate is here
        uint64_t prime_candidate = at_least;
        // Check if the given number is less than 2
        if (prime_candidate <= 2)
            return 2;
        // If the given number is even it cannot be a prime so we make it odd
        if (prime_candidate % 2 == 0)
            prime_candidate += 1;
        // Define few useful variables
        uint64_t max_check;
        uint64_t check;
        bool is_prime;
        // Run this loop until a prime is found
        while(true)
        {
            // Prime divisibility checking starts at 3
            check = 3;
            // Prime divisibility checking ends at the square root
            max_check = std::floor(std::sqrt(prime_candidate));
            is_prime=true;
            //std::cout << "Checking for prime "<< prime_candidate << " with max check " << max_check << "\n";
            // Check divisibility for all relevant values
            while(check <= max_check)
            {
                // If divisible, no prime
                if (prime_candidate % check == 0)
                {
                    is_prime = false;
                    break;
                }
                check+=1;
            }
            // When we find a prime, rreturn it
            if (is_prime)
            {
                std::cout << "Hash table size is: " << prime_candidate << "\n";
                return prime_candidate;
            }
                
            // If candidate is not prime, increase it by 2 (skip the even value)
            prime_candidate+=2;
        }
    }

    uint64_t next_prime3mod4(uint64_t at_least)
    {
        // The prime candidate is here
        uint64_t prime_candidate = at_least;
        // Check if the given number is less than 2
        if (prime_candidate <= 2)
            return 2;
        // If the given number is even it cannot be a prime so we make it odd
        if (prime_candidate % 2 == 0)
            prime_candidate += 1;
        // Define few useful variables
        uint64_t max_check;
        uint64_t check;
        bool is_prime;
        // Run this loop until a prime is found
        while(true)
        {
            // Prime divisibility checking starts at 3
            check = 3;
            // Prime divisibility checking ends at the square root
            max_check = std::floor(std::sqrt(prime_candidate));
            is_prime=true;
            //std::cout << "Checking for prime "<< prime_candidate << " with max check " << max_check << "\n";
            // Check divisibility for all relevant values
            while(check <= max_check)
            {
                // If divisible, no prime
                if (prime_candidate % check == 0)
                {
                    is_prime = false;
                    break;
                }
                check+=1;
            }
            // When we find a prime, rreturn it
            if ((is_prime) && (prime_candidate % 4 == 3))
            {
                std::cout << "Hash table size is: " << prime_candidate << "\n";
                return prime_candidate;
            }
            // If candidate is not prime, increase it by 2 (skip the even value)
            prime_candidate+=2;
        }
    }


    uint64_t modular_multiplicative_inverse_coprimes(int64_t A, int64_t M)
    {
        //int64_t AA = A;
        //int64_t MM = M;
        int64_t m0 = M;
        int64_t y = 0, x = 1;
    
        if (M == 1)
            return 0;
    
        while (A > 1) {
            // q is quotient
            int64_t q = A / M;
            int64_t t = M;
    
            // m is remainder now, process same as
            // Euclid's algo
            M = A % M, A = t;
            t = y;
    
            // Update y and x
            y = x - q * y;
            x = t;
        }
    
        // Make x positive
        if (x < 0)
            x += m0;
    
        //std::cout << "Modular multiplicative inverse of A=" << AA << " and M=" << MM << " is " << x << "\n";
        return uint64_t(x);
    }

    uint64_t modular_multiplicative_inverse_little_fermat(uint64_t A, uint64_t M)
    {
        /*
        bool nbf = false;
        for (int X = 1; X < M; X++){
            if (((A % M) * (X % M)) % M == 1){
                std::cout << "Naive brute force modular multiplicative inverse of A=" << A << " and M=" << M << " is " << X << "\n";
                nbf=true;
                break;
                //return X;
            }
        }
        if (!nbf){
            std::cout << "Modular multiplicative inverse for A=" << A << " and M=" << M << " does not exist\n";
        }
        */
            
        uint64_t g = greatest_common_divisor(A, M);
        if (g != 1)
        {
            std::cout << "Modular multiplicative inverse does not exist. Invalid hash table size M.";
            exit(1);
        }
        return power_under_modulo(A, M - 2, M);
    }

    uint64_t power_under_modulo(uint64_t x, uint64_t y, uint64_t M)
    {
        if (y == 0)
            return 1;
        uint64_t p = power_under_modulo(x, y/2, M) % M;
        p = (p * p) % M;

        return (y % 2 == 0) ? p : (x * p) % M;
    }

    // Smaller value first?
    uint64_t greatest_common_divisor(uint64_t A, uint64_t B)
    {
        if (A == 0)
            return B;
        return greatest_common_divisor(B % A, A);
    }
}