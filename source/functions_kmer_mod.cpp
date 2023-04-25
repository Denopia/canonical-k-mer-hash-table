#include "functions_kmer_mod.hpp"

namespace kmod
{
    //===================
    // Modifier functions
    //===================

    uint64_t modify_to_increase_count_by_one(uint64_t D)
    {
        uint64_t D2 = D;
        if (((D2 >> 12) & uint64_t(16383)) != uint64_t(16383))
            D2 = D2 + uint64_t(4096);
        else
            std::cout << "Count was not increased\n";
        return D2;
    }


    uint64_t modify_predecessor_slot(uint64_t D, uint64_t predecessor_slot)
    {
        uint64_t D2 = D;
        D2 = (D2 & uint64_t(67108863));
        D2 = (D2 | (predecessor_slot<<26));
        return D2;
    }

    uint64_t modify_left_character(uint64_t D, uint64_t left_char)
    {
        uint64_t D2 = D;
        D2 = (D2 & (~(uint64_t(3072))));
        D2 = (D2 | (left_char << 10));
        return D2;
    }

    uint64_t modify_right_character(uint64_t D, uint64_t right_char)
    {
        uint64_t D2 = D;
        D2 = (D2 & (~(uint64_t(768))));
        D2 = (D2 | (right_char << 8));
        return D2;
    }

    uint64_t modify_to_occupied(uint64_t D)
    {
        uint64_t D2 = D;
        D2 = D2 | (uint64_t(1));
        return D2;
    }

    uint64_t modify_to_unoccupied(uint64_t D)
    {
        uint64_t D2 = D;
        D2 = D2 & ~(uint64_t(1));
        return D2;
    }

    uint64_t modify_to_have_predecessor(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<1);
        return D2;
    }

    uint64_t modify_to_not_have_predecessor(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 1);
        return D2;
    }

    uint64_t modify_to_left_char_is_null(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<2);
        return D2;
    }

    uint64_t modify_to_left_char_is_notnull(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 2);
        return D2;
    }

    uint64_t modify_to_right_char_is_null(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<3);
        return D2;
    }

    uint64_t modify_to_right_char_is_notnull(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 3);
        return D2;
    }

    uint64_t modify_to_be_canonical_during_insertion_self(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<4);
        return D2;
    }

    uint64_t modify_to_be_noncanonical_during_insertion_self(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 4);
        return D2;
    }

    uint64_t modify_to_be_canonical_during_insertion_pred(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<5);
        return D2;
    }

    uint64_t modify_to_be_noncanonical_during_insertion_pred(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 5);
        return D2;
    }

    uint64_t modify_to_be_flagged_1(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<6);
        return D2;
    }

    uint64_t modify_to_be_unflagged_1(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 6);
        return D2;
    }

    uint64_t modify_to_be_flagged_2(uint64_t D)
    {
        uint64_t D2 = D;
        D2 |= (uint64_t(1)<<7);
        return D2;
    }

    uint64_t modify_to_be_unflagged_2(uint64_t D)
    {
        uint64_t D2 = D;
        D2 &= ~(uint64_t(1) << 7);
        return D2;
    }

    uint64_t modify_for_migration(uint64_t D, uint64_t left_char, uint64_t right_char, uint64_t predecessor_slot, uint64_t self_canonical_during_insertion, uint64_t pred_canonical_during_insertion)
    {
        uint64_t D2 = D;
        D2 = kmod::modify_left_character(D2, left_char);
        D2 = kmod::modify_right_character(D2, right_char);
        D2 = kmod::modify_to_have_predecessor(D2);
        D2 = kmod::modify_predecessor_slot(D2, predecessor_slot);
        if (self_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_self(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_self(D2);
        }
        if (pred_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_pred(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_pred(D2);
        }
        return D2;   
    }

    uint64_t modify_for_insertion(uint64_t D, bool predecessor_exists, bool pred_canonical_during_insertion, uint64_t predecessor_slot, bool self_canonical_during_insertion, uint64_t left_char, uint64_t right_char)
    {
        uint64_t D2 = 0;
        D2 = kmod::modify_to_occupied(D2);
        D2 = kmod::modify_to_increase_count_by_one(D2);
        D2 = kmod::modify_left_character(D2, left_char);
        D2 = kmod::modify_right_character(D2, right_char);
        if (predecessor_exists){
            D2 = kmod::modify_to_have_predecessor(D2);
        } else {
            D2 = kmod::modify_to_not_have_predecessor(D2);
        }
        D2 = kmod::modify_predecessor_slot(D2, predecessor_slot);
        if (self_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_self(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_self(D2);
        }
        if (pred_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_pred(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_pred(D2);
        }
        return D2;   
    }

    uint64_t modify_predecessor_slot_and_orientations(uint64_t D, uint64_t predecessor_slot, bool self_canonical_during_insertion, bool pred_canonical_during_insertion)
    {
        uint64_t D2 = D;
        D2 = kmod::modify_predecessor_slot(D2, predecessor_slot);
        if (self_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_self(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_self(D2);
        }
        if (pred_canonical_during_insertion){
            D2 = kmod::modify_to_be_canonical_during_insertion_pred(D2);
        } else {
            D2 = kmod::modify_to_be_noncanonical_during_insertion_pred(D2);
        }
        return D2;
    }

}