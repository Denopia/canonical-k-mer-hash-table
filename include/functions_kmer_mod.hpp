#include <cstdint>
#include <iostream>

#pragma once

namespace kmod
{
    // Data modifiers, given uint64_t D and the desired operation, return D after operation
    uint64_t modify_to_increase_count_by_one(uint64_t D);
    uint64_t modify_predecessor_slot(uint64_t D, uint64_t predecessor_slot);
    uint64_t modify_predecessor_slot_and_orientations(uint64_t D, uint64_t predecessor_slot, bool self_canonical_during_insertion, bool pred_canonical_during_insertion);
    uint64_t modify_left_character(uint64_t D, uint64_t left_char);
    uint64_t modify_right_character(uint64_t D, uint64_t right_char);
    uint64_t modify_to_occupied(uint64_t D);
    uint64_t modify_to_unoccupied(uint64_t D);
    uint64_t modify_to_have_predecessor(uint64_t D);
    uint64_t modify_to_not_have_predecessor(uint64_t D);
    uint64_t modify_to_left_char_is_null(uint64_t D);
    uint64_t modify_to_left_char_is_notnull(uint64_t D);
    uint64_t modify_to_right_char_is_null(uint64_t D);
    uint64_t modify_to_right_char_is_notnull(uint64_t D);
    uint64_t modify_to_be_canonical_during_insertion_self(uint64_t D);
    uint64_t modify_to_be_noncanonical_during_insertion_self(uint64_t D);
    uint64_t modify_to_be_canonical_during_insertion_pred(uint64_t D);
    uint64_t modify_to_be_noncanonical_during_insertion_pred(uint64_t D);
    uint64_t modify_to_be_flagged_1(uint64_t D);
    uint64_t modify_to_be_unflagged_1(uint64_t D);
    uint64_t modify_to_be_flagged_2(uint64_t D);
    uint64_t modify_to_be_unflagged_2(uint64_t D);
    uint64_t modify_for_migration(uint64_t D, uint64_t left_char, uint64_t right_char, uint64_t  predecessor_slot, uint64_t self_forward_canonical, uint64_t pred_forward_canonical);
    uint64_t modify_for_insertion(uint64_t D, bool predecessor_exists, bool pred_canonical_during_insertion, uint64_t predecessor_slot, bool self_canonical_during_insertion, uint64_t left_char, uint64_t right_char);
}