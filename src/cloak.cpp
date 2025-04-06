/*
 * CLOAK
 * Originally authored by Chiragdeep Chatur (Masel Lab, University of Arizona)
 * Translated into C++ and integrated into MUSCLE by Peter Goodman (Masel Lab, University of Arizona)
 * Description: CLOAK takes in a set of variant alignments of the same sequences and searches for 
 * consensus among them. 
 * To generate the required input alignments for CLOAK, use Muscle5 options to conduct three 
 * perturbations of its core Hidden Markov Model (HMM) and three perturbations of its core 
 * guide tree, the combination of which produces 16 alternative alignments. 
 * The CLOAK algorithm retains each pair of amino acids in the same column if and only if they are 
 * found together (in the same column) in all 16 alternative alignments. 
 * If multiple subsets of amino acids within a column are found together across all alignments, 
 * but the whole column is not, the subsets are split into multiple columns that are all retained.
 *
 * Usage:       muscle -cloak input_ensemble_file -mincol <integer> -output <output_file_name>
 *
 * Arguments:
 *    - input_ensemble_file : Path to the input MSA file, which can either be an EFA file 
 *                            or a text file with paths to individual MSAs on each line
 *    - -mincol <integer>   : Minimum number of non-gap characters required per column
 *                            for that column to be retained in the output.
 *                            Default value of 2 if not specified
 *    - -output <filename>  : Name of the file where the filtered MSA will be written.
 *
 * Notes:
 *    - This tool is intended to be used as a plug-in with the MUSCLE suite.
 *    - Columns not meeting the minimum occupancy threshold will be discarded.
 *    - Input ensemble must either be EFA file or text file with paths to individual MSAs
 *
 * Example:
 *    muscle -cloak input_aligmments.efa -mincol 5 -output cleaned_alignment.fasta
 *
 */

#include "muscle.h"
#include "myutils.h"
#include "ensemble.h"
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <cassert>
#include <fstream>


//Converts the users input MSA ensemble into a 3-D array
//Takes either an EFA of text file
//3-D array format: MSA->Row->Char
std::vector<std::vector<std::vector<char>>> convertEnsembleTo3D(
    const Ensemble& E,
    std::vector<std::vector<std::string>>& allLabels)
{
    // Build a 3D array: each MSA => 2D array => rows of chars
    std::vector<std::vector<std::vector<char>>> allAlignments;

    uint nMSAs = E.GetMSACount();
    allLabels.resize(nMSAs);

    for (uint i = 0; i < nMSAs; i++)
    {
        const MSA& m = E.GetMSA(i);
        uint seqCount = m.GetSeqCount();
        uint colCount = m.GetColCount();

        // Gather (label, sequence) pairs, then sort by label
        struct LabelSeq
        {
            std::string label;
            std::vector<char> seq;
        };
        std::vector<LabelSeq> rowData;
        rowData.reserve(seqCount);

        for (uint s = 0; s < seqCount; s++)
        {
            LabelSeq ls;
            ls.label = m.GetLabel(s);
            ls.seq.reserve(colCount);

            for (uint c = 0; c < colCount; c++)
            {
                ls.seq.push_back(m.GetChar(s, c));
            }
            rowData.push_back(ls);
        }

        // Sort by label
        std::sort(rowData.begin(), rowData.end(),
            [](const LabelSeq& a, const LabelSeq& b) {
                return a.label < b.label;
            });

        // Build 2D char array
        std::vector<std::vector<char>> twoD;
        twoD.reserve(seqCount);

        std::vector<std::string> labelVec;
        labelVec.reserve(seqCount);

        for (auto& ls : rowData)
        {
            twoD.push_back(ls.seq);
            labelVec.push_back(ls.label);
        }

        allAlignments.push_back(twoD);
        allLabels[i] = labelVec;
    }

    return allAlignments; // [MSA_index][row_index][char_index]
}

// Removes leading and trailing whitespace from a string.
std::string trim(const std::string& s) {
    std::string result = s;
    result.erase(result.begin(),
        std::find_if(result.begin(), result.end(),
            [](unsigned char ch) { return !std::isspace(ch); }));
    result.erase(std::find_if(result.rbegin(), result.rend(),
        [](unsigned char ch) { return !std::isspace(ch); }).base(),
        result.end());
    return result;
}


// Helper class/structure/functions to deal with mixed vectors
// We need to convert sequence vectors from letters to ints
// Dashes must stay as Chars, so we use a helper class
enum class NumOrDashType { INTEGER, CHARACTER };

struct NumOrDash {
    NumOrDashType type;
    union {
        int intValue;
        char charValue;
    };

    // Actual constructor for int
    NumOrDash(int v) {
        type = NumOrDashType::INTEGER;
        intValue = v;
    }

    // Actual constructor for char
    NumOrDash(char c) {
        type = NumOrDashType::CHARACTER;
        charValue = c;
    }
};

//   For each alignment (a 2-D vector of characters) in the 3-D array,
//   and for each sequence (row) in that alignment, assign a running count
//   (starting at 1) to each non-dash character, leaving dashes unchanged.
//   The result is returned as a 3-D vector, where each element is a NumOrDash.
std::vector<std::vector<std::vector<NumOrDash>>>
convert_to_nums(const std::vector<std::vector<std::vector<char>>>& array)
{
    std::vector<std::vector<std::vector<NumOrDash>>> nums;

    // Iterate over each alignment (each 2-D array)
    for (size_t x = 0; x < array.size(); ++x) {
        const auto& two_d_array = array[x];
        std::vector<std::vector<NumOrDash>> make_two_d_array;

        // Process each sequence (row) in the alignment.
        for (size_t y = 0; y < two_d_array.size(); ++y) {
            const auto& row = two_d_array[y];
            int i = 1;  // Reset counter for each sequence.
            std::vector<NumOrDash> row_array;

            // Process each character in the sequence.
            for (size_t z = 0; z < row.size(); ++z) {
                char letter = row[z];
                if (letter == '-') {
                    // Calls the NumOrDash(char) constructor
                    row_array.push_back(NumOrDash('-'));
                }
                else {
                    // Calls the NumOrDash(int) constructor
                    row_array.push_back(NumOrDash(i));
                    i++;
                }
            }
            make_two_d_array.push_back(row_array);
        }
        nums.push_back(make_two_d_array);
    }
    return nums;
}


// toSearch: A 1-D vector of NumOrDash representing the column to match.
// three_d_array_of_divvies: A 3-D vector, where each element is an alignment
//   (a 2-D vector of rows) and each row is a vector of NumOrDash.
// current_alignment_number: The index of the alignment that should be skipped.
// For each alignment (except the one specified by current_alignment_number),
// the function checks whether there is at least one row (divvy) such that for every
// index k, if toSearch[k] is not a dash then toSearch[k] equals divvy[k].
// Returns true if every alignment (other than the current one) has at least one matching row
// otherwise, false.
bool can_this_divvy_be_found_without_subsequent_divvies(
    const std::vector<NumOrDash>& toSearch,
    const std::vector<std::vector<std::vector<NumOrDash>>>& three_d_array_of_divvies,
    size_t current_alignment_number)
{
    // Iterate over each alignment.
    for (size_t x = 0; x < three_d_array_of_divvies.size(); ++x) {
        if (x == current_alignment_number)
            continue;

        const auto& divvies_of_current_alignment = three_d_array_of_divvies[x];
        bool current_alignment_has_match = false;

        // For each "row" in this alignment
        for (const auto& divvy : divvies_of_current_alignment) {
            bool match_was_found = true;

            // Compare each element in the toSearch vector.
            for (size_t k = 0; k < toSearch.size(); ++k) {
                const NumOrDash& current_num = toSearch[k];
                const NumOrDash& compare_to = divvy[k];

                // If current_num is a dash, skip the comparison.
                // (type == CHARACTER && charValue == '-')
                if (current_num.type == NumOrDashType::CHARACTER &&
                    current_num.charValue == '-')
                {
                    // We do nothing (skip).
                    continue;
                }
                else {
                    // Otherwise, first ensure the types match.
                    if (current_num.type != compare_to.type) {
                        match_was_found = false;
                        break;
                    }

                    // Compare the values depending on type.
                    if (current_num.type == NumOrDashType::INTEGER) {
                        if (current_num.intValue != compare_to.intValue) {
                            match_was_found = false;
                            break;
                        }
                    }
                    else { // NumOrDashType::CHARACTER
                        if (current_num.charValue != compare_to.charValue) {
                            match_was_found = false;
                            break;
                        }
                    }
                }
            } // end for each element k

            if (match_was_found) {
                current_alignment_has_match = true;
                break; // No need to check more rows in this alignment
            }
        } // end for each divvy

        // If no row in this alignment matched, return false
        if (!current_alignment_has_match) {
            return false;
        }
    }

    return true;
}


// For each dictionary in the input array (each representing one alignment),
// iterate over each key (each with an associated vector of NumOrDash).
// Count the non-dash elements in that vector, and if there are at least 2,
// include that vector in the alignment's output. Finally, return a 3-D array
// (vector of alignments, where each alignment is a vector of rows, and each row
// is a vector of NumOrDash).                  
std::vector<std::vector<std::vector<NumOrDash>>>
divvies_as_3darray(const std::vector<std::map<int, std::vector<NumOrDash>>>& array_of_dicts)
{
    std::vector<std::vector<std::vector<NumOrDash>>> return_array;

    // Iterate over each dictionary (each alignment).
    for (const auto& divvy_dict : array_of_dicts) {
        std::vector<std::vector<NumOrDash>> array_per_alignment;

        // For each key-value pair in the dictionary:
        for (const auto& entry : divvy_dict) {
            const auto& curr_array = entry.second;
            int number_of_nums = 0;

            // Count the number of elements that are not dashes.
            for (const auto& elem : curr_array) {
                // If elem is an integer, it's never a dash.
                if (elem.type == NumOrDashType::INTEGER) {
                    number_of_nums++;
                }
                else {
                    // If it's a character, check if it's not '-'.
                    if (elem.charValue != '-') {
                        number_of_nums++;
                    }
                }
            }

            // If at least 2 non-dash elements exist, include this vector.
            if (number_of_nums >= 2) {
                array_per_alignment.push_back(curr_array);
            }
        }

        return_array.push_back(array_per_alignment);
    }

    return return_array;
}

// For each indices array in all_indices, build a dictionary (map) where:
//   - Each key is a non-dash number found in the indices array.
//   - Then for each position in the indices array:
//       * If the element is a dash, append '-' to every key�s vector.
//       * Otherwise, let curr be the number; then append to myDict[curr] the
//         corresponding element from toSearch, and for every other key append '-'.
// Returns a vector of such dictionaries.
std::vector<std::map<int, std::vector<NumOrDash>>>
get_divvied_results_as_array_of_dicts(const std::vector<std::vector<NumOrDash>>& all_indices,
    const std::vector<NumOrDash>& toSearch)
{
    std::vector<std::map<int, std::vector<NumOrDash>>> someArray;

    // Process each "indices" array.
    for (const auto& indices : all_indices) {
        std::map<int, std::vector<NumOrDash>> myDict;

        // First pass: add a key for every non-dash element.
        for (size_t i = 0; i < indices.size(); i++) {
            // If indices[i] holds an int, that int is our key.
            if (indices[i].type == NumOrDashType::INTEGER) {
                int key = indices[i].intValue;
                if (myDict.find(key) == myDict.end()) {
                    myDict[key] = std::vector<NumOrDash>();
                }
            }
        }

        // Second pass: fill in the vectors.
        for (size_t i = 0; i < indices.size(); i++) {
            // If the element is a dash (char == '-'), append '-' to every key's vector.
            if (indices[i].type == NumOrDashType::CHARACTER &&
                indices[i].charValue == '-')
            {
                for (auto& pair : myDict) {
                    pair.second.push_back(NumOrDash('-'));
                }
            }
            else {
                // Otherwise, the element is an integer key.
                int curr = indices[i].intValue;

                // Append the corresponding toSearch element for the matching key.
                myDict[curr].push_back(toSearch[i]);

                // For every other key, append '-'.
                for (auto& pair : myDict) {
                    if (pair.first != curr) {
                        pair.second.push_back(NumOrDash('-'));
                    }
                }
            }
        }

        someArray.push_back(myDict);
    }

    return someArray;
}

// For each indices array in 'all_indices', this function skips any leading
// dashes, then takes the first non-dash value as a reference number. It then
// examines the remainder of the array (ignoring dashes). If any non-dash value
// is different from the reference, it returns true (divvy is required).
// If no indices array shows a discrepancy, it returns false.
bool divvy_is_required(const std::vector<std::vector<NumOrDash>>& all_indices)
{
    for (const auto& indices_array : all_indices) {
        if (indices_array.empty()) {
            continue;
        }

        size_t i = 0;
        // Skip over any leading dashes.
        while (i < indices_array.size() &&
            indices_array[i].type == NumOrDashType::CHARACTER &&
            indices_array[i].charValue == '-')
        {
            i++;
        }

        // If the entire array is dashes, move on.
        if (i == indices_array.size()) {
            continue;
        }

        // The first non-dash element should be an int.
        int number_to_match_with = 0;
        if (indices_array[i].type == NumOrDashType::INTEGER) {
            number_to_match_with = indices_array[i].intValue;
        }
        else {
            // Data unexpected, skip to next.
            continue;
        }

        // Check the remainder of indices_array.
        for (size_t x = i; x < indices_array.size(); ++x) {
            // If the element is a dash, skip it.
            if (indices_array[x].type == NumOrDashType::CHARACTER &&
                indices_array[x].charValue == '-')
            {
                continue;
            }
            // Else if the element is an integer, compare.
            else if (indices_array[x].type == NumOrDashType::INTEGER) {
                int currentNum = indices_array[x].intValue;
                if (currentNum != number_to_match_with) {
                    return true; // Mismatch found
                }
            }
        }
    }
    return false; // No mismatch found anywhere
}


// Given the 3-D array of alignments (all_alignments), a 1-D array (toSearch)
// containing one element from each sequence of the first alignment, and an integer
// (curr_alignment_number) indicating which alignment to process, this function
// searches each row of the specified alignment for the element from toSearch.
// If the element is a dash ('-'), a dash is appended to the result;
// otherwise, for each occurrence in the row equal to the element, the index is appended.
std::vector<NumOrDash> create_indices(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_alignments,
    const std::vector<NumOrDash>& toSearch,
    size_t curr_alignment_number)
{
    std::vector<NumOrDash> indices;
    const auto& current_alignment = all_alignments[curr_alignment_number];

    for (size_t x = 0; x < current_alignment.size(); ++x) {
        const auto& row = current_alignment[x];
        NumOrDash num = toSearch[x];

        // If the toSearch element is a dash, record a dash.
        if (num.type == NumOrDashType::CHARACTER && num.charValue == '-') {
            indices.push_back(NumOrDash('-'));
        }
        else {
            // Otherwise, for each element in the row, if it equals num, append its index.
            for (size_t a = 0; a < row.size(); ++a) {
                // Compare only if the types match.
                if (row[a].type == num.type) {
                    if (row[a].type == NumOrDashType::INTEGER) {
                        if (row[a].intValue == num.intValue) {
                            // Push back 'a' as a NumOrDash integer
                            indices.push_back(NumOrDash(static_cast<int>(a)));
                        }
                    }
                    else { // NumOrDashType::CHARACTER
                        if (row[a].charValue == num.charValue) {
                            indices.push_back(NumOrDash(static_cast<int>(a)));
                        }
                    }
                }
            }
        }
    }
    return indices;
}


// Helper function: creates an initial 1D vector to search
std::vector<NumOrDash> create_toSearch(const std::vector<std::vector<NumOrDash>>& alignment, size_t col) {
    std::vector<NumOrDash> result;
    for (const auto& row : alignment) {
        if (col < row.size())
            result.push_back(row[col]);
    }
    return result;
}

// Helper function, takes in a 2D Array, transposes it
template <typename Container>
std::vector<std::vector<typename Container::value_type>> transpose_output(const std::vector<Container>& l1) {
    std::vector<std::vector<typename Container::value_type>> l2;
    if (l1.empty()) return l2;
    size_t cols = l1[0].size();
    for (size_t i = 0; i < cols; ++i) {
        std::vector<typename Container::value_type> row;
        for (const auto& item : l1) {
            if (i < item.size()) {
                row.push_back(item[i]);
            }
            else {
                break;
            }
        }
        l2.push_back(row);
    }
    return l2;
}


// This function takes a vector of sequences (strings) and returns a new vector 
// containing only those sequences that have more than one non-dash character.
std::vector<std::vector<char>> remove_singletons(
    const std::vector<std::vector<char>>& two_d_array,
    uint minColumns)
{
    std::vector<std::vector<char>> return_array;
    for (const auto& one_d_array : two_d_array) {
        unsigned int num_letters = 0;
        for (char letter : one_d_array) {
            if (letter != '-') {
                ++num_letters;
            }
        }
        if (num_letters >= minColumns) {
            return_array.push_back(one_d_array);
        }
    }
    return return_array;
}


// Checks whether two vectors of NumOrDash objects are equal
// in both size and content.Two NumOrDash elements are considered equal if:
// Their types(NumOrDashType::INTEGER or CHARACTER) match.
// Their corresponding values(intValue or charValue) are equal.
// The function is used, for example, to eliminate duplicate NumOrDash vectors
// from a collection.

bool equalNumOrDashVector(const std::vector<NumOrDash>& a,
    const std::vector<NumOrDash>& b)
{
    if (a.size() != b.size()) {
        return false;
    }

    for (size_t i = 0; i < a.size(); i++) {
        // Check if types match
        if (a[i].type != b[i].type) {
            return false;
        }

        // Compare appropriate union member
        if (a[i].type == NumOrDashType::INTEGER) {
            if (a[i].intValue != b[i].intValue) {
                return false;
            }
        }
        else { // CHARACTER
            if (a[i].charValue != b[i].charValue) {
                return false;
            }
        }
    }
    return true;
}


// Recursively filters and adds valid divvied alignments to the output set.
// Given an initial vector of alignment indices(`toSearch`) from the first alignment,
// this function explores all alternative "divvied" representations derived from the
// other alignments and determines whether each should be added to the output.
// A candidate divvy is added to the output if it can be found without relying on
// subsequent divvies(i.e., it is valid independently).Otherwise, the function
// recursively evaluates the dependencies of that divvy before deciding.
// Divvies identical to the original `toSearch` are skipped to avoid redundant work.

void add_to_output_or_discard(
    const std::vector<std::vector<NumOrDash>>& indices_of_all_alignments,
    std::vector<std::vector<NumOrDash>>& output,
    const std::vector<NumOrDash>& toSearch,
    std::vector<std::vector<NumOrDash>>& seen)
{
    // Count the number of non-dash elements in toSearch.
    int number_of_nums = 0;
    for (const auto& num : toSearch) {
        if (num.type == NumOrDashType::INTEGER) {
            number_of_nums++;
        }
        else {
            // It's a character; check if it's not '-'
            if (num.charValue != '-') {
                number_of_nums++;
            }
        }
    }

    if (number_of_nums == 0) {
        return;
    }

    auto array_containing_dicts_of_divvies =
        get_divvied_results_as_array_of_dicts(indices_of_all_alignments, toSearch);
    auto three_d_array_of_divvies =
        divvies_as_3darray(array_containing_dicts_of_divvies);

    // Loop starting from index 1 (skip alignment 0).
    for (size_t i = 1; i < three_d_array_of_divvies.size(); i++) {
        const auto& divvies_from_specific_alignment = three_d_array_of_divvies[i];
        for (const auto& divvy : divvies_from_specific_alignment) {
            // If toSearch equals divvy, skip it.
            // (Assumes you have a function or operator== that compares vectors of NumOrDash)
            if (equalNumOrDashVector(divvy, toSearch)) {
                continue;
            }
            bool can_found = can_this_divvy_be_found_without_subsequent_divvies(
                divvy, three_d_array_of_divvies, i);

            if (can_found) {
                output.push_back(divvy);
            }
            else {
                add_to_output_or_discard(indices_of_all_alignments, output, divvy, seen);
            }
        }
    }
}


// Constructs the initial numeric consensus output from aligned input sequences.
// For each column in the first alignment of the ensemble, this function generates a
// representative column(`toSearch`) and checks whether it is consistent across all
// alignments.If consistent, it is directly added to the output.If not, the column
// is "divvied" into one or more alternative representations based on alignment - specific
// mappings.Valid divvied representations are recursively added via
// add_to_output_or_discard().
//
// param all_als A 3D vector representing all input alignments in NumOrDash format,
// where all_als[i][r][c] corresponds to column c in row r of alignment i.
// return A 2D vector of NumOrDash elements representing the consensus output.
std::vector<std::vector<NumOrDash>> create_output(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_als)
{
    std::vector<std::vector<NumOrDash>> output;
    const auto& first_alignment = all_als[0];
    // Use the number of elements in the first row as the number of columns.
    size_t num_cols = first_alignment[0].size();

    for (size_t x = 0; x < num_cols; x++) {
        // Create the initial toSearch term (i.e. column x of first_alignment)
        std::vector<NumOrDash> toSearch = create_toSearch(first_alignment, x);

        // Build indices_of_all_alignments: one indices array per alignment.
        std::vector<std::vector<NumOrDash>> indices_of_all_alignments;
        for (size_t k = 0; k < all_als.size(); k++) {
            std::vector<NumOrDash> indices = create_indices(all_als, toSearch, k);
            indices_of_all_alignments.push_back(indices);
        }

        // Check if divvy is needed.
        bool divvy_needed = divvy_is_required(indices_of_all_alignments);
        if (!divvy_needed) {
            output.push_back(toSearch);
        }
        else {
            std::vector<std::vector<NumOrDash>> seen; // initially empty
            add_to_output_or_discard(indices_of_all_alignments, output, toSearch, seen);
        }
        // Reset toSearch is not necessary.
    }
    return output;
}

// nums_array: a 2D vector (each inner vector is a row) of NumOrDash (each element is either an int or '-' as a char)
// first_alignment: a 2D vector of char (all letters)
// Returns a 2D vector of char.
std::vector<std::vector<char>> produce_output_in_letters(
    const std::vector<std::vector<NumOrDash>>& nums_array,
    const std::vector<std::vector<char>>& first_alignment)
{
    std::vector<std::vector<char>> return_array;
    std::vector<std::vector<char>> intermediary;

    // If first_alignment is empty, there's nothing to do
    if (first_alignment.empty() || first_alignment[0].empty()) {
        return return_array;
    }

    // Build an "intermediary" 2D array that collects the non-dash
    // characters column by column from first_alignment.
    size_t num_cols = first_alignment[0].size();
    for (size_t i = 0; i < num_cols; i++) {
        std::string s;
        // Concatenate the i-th element from every row
        for (const auto& row : first_alignment) {
            if (i < row.size()) {
                s.push_back(row[i]);
            }
        }
        // Build an array of characters from s, skipping dashes.
        std::vector<char> some_array;
        for (char ch : s) {
            if (ch == '-') {
                continue;
            }
            some_array.push_back(ch);
        }
        intermediary.push_back(some_array);
    }

    // For each row in nums_array:
    for (const auto& nums : nums_array) {
        std::vector<char> subarray;
        size_t cols = nums.size();

        for (size_t i = 0; i < cols; i++) {
            // If it is a dash, we output '-'.
            if (nums[i].type == NumOrDashType::CHARACTER &&
                nums[i].charValue == '-')
            {
                subarray.push_back('-');
            }
            // Else if it's an integer:
            else if (nums[i].type == NumOrDashType::INTEGER) {
                int int_val = nums[i].intValue;
                int index = int_val - 1;  // old code: int index = std::get<int>(nums[i]) - 1;

                // Check bounds in intermediary[i]
                if (i < intermediary.size() && index >= 0 &&
                    index < static_cast<int>(intermediary[i].size()))
                {
                    subarray.push_back(intermediary[i][index]);
                }
                else {
                    // If out-of-range, append a dash.
                    subarray.push_back('-');
                }
            }
            else {
                // Fallback case
                subarray.push_back('-');
            }
        }
        return_array.push_back(subarray);
    }
    return return_array;
}


void cmd_cloak()

{
    // Parse user options
    std::string ensemblePath = opt(cloak);  // user runs: muscle -cloak myfile ...
    if (ensemblePath.empty())
    {
        Log("Error: you must specify -cloak <ensembleFile or list-of-FASTA>\n");
        return;
    }

    std::string outputFilename = opt(output);
    if (outputFilename.empty())
        outputFilename = "cloak_result.fasta";

    uint minimumColumns = 2; // default if -mincol is not specified
    uint userMincol = opt(mincol);  // Parser returns 0 if not specified

    // If the parser's sentinel is 0 when -mincol isn't used, then check:
    if (userMincol != 0)
    {
        // User actually specified -mincol
        minimumColumns = userMincol;
    }
    if (minimumColumns < 2)
    {
        Log("Error: -mincol must be >= 2\n");
        return;
    }

    Log("Starting cloak. Ensemble path: %s, Output filename: %s, Minimum columns requested: %u\n",
        ensemblePath.c_str(), outputFilename.c_str(), userMincol);

    // Load the entire ensemble from the user-supplied file
    Ensemble E;
    try
    {
        E.FromFile(ensemblePath);
    }
    catch (...)
    {
        Log("Error: Unable to parse ensemble from file %s\n", ensemblePath.c_str());
        return;
    }

    // Convert the ensemble into a 3D array of chars
    // plus keep a parallel 2D array of labels for each MSA.
    // allLabels[i][r] = label for row r in the i-th MSA

    Log("Converting ensemble to 3D array of letters...\n");

    std::vector<std::vector<std::string>> allLabels;
    auto three_d_array_letters = convertEnsembleTo3D(E, allLabels);


    // CLOAK Algorithm
    // 1) Convert letters -> numbers (and dashes)
    Log("Converting letters to numeric representation.\n");
    auto nums_answer = convert_to_nums(three_d_array_letters);

    // 2) Create initial numeric output
    Log("Creating initial numeric output.\n");
    auto nums_array_output_initial = create_output(nums_answer);

    // 3) Remove duplicates
    Log("Removing duplicates from numeric output.\n");
    std::vector<std::vector<NumOrDash>> nums_array_output;
    for (const auto& num_array : nums_array_output_initial)
    {
        bool duplicate = false;
        for (const auto& existing : nums_array_output)
        {
            if (equalNumOrDashVector(existing, num_array))
            {
                duplicate = true;
                break;
            }
        }
        if (!duplicate)
            nums_array_output.push_back(num_array);
    }

    // 4) Use the first alignment in three_d_array_letters as reference.
    auto reference_alignment = three_d_array_letters[0];

    // 5) Convert numeric output back to letters
    Log("Converting numeric output back to letters.\n");
    auto transposed_reference = transpose_output(reference_alignment);
    auto letters_array_output_initial =
        produce_output_in_letters(nums_array_output, transposed_reference);

    Log("Removing singletons with minimum column threshold = %u.\n", minimumColumns);
    auto letters_array_output = remove_singletons(letters_array_output_initial, minimumColumns);

    letters_array_output = transpose_output(letters_array_output);

    // Build final MSA with original labels from the original MSAs

    // allLabels[0] = the labels (in alphabetical order) for the first MSA
    const auto& firstLabels = allLabels[0];

    // Build finalLabels/finalSeqs in parallel
    std::vector<std::string> finalLabels;
    std::vector<std::string> finalSeqs;

    for (size_t i = 0; i < letters_array_output.size(); i++)
    {
        // Convert the row to a string
        std::string rowStr(letters_array_output[i].begin(),
            letters_array_output[i].end());
        // Count non-dash letters
        int letterCount = 0;
        for (char c : rowStr)
            if (c != '-') letterCount++;
        // If this row is NOT a singleton, we keep it
        if (letterCount > 1)
        {
            // The label for row i is firstLabels[i]
            if (i < firstLabels.size())
                finalLabels.push_back(firstLabels[i]);
            else
                finalLabels.push_back(std::string(">Seq") + std::to_string(i + 1));

            finalSeqs.push_back(rowStr);
        }
    }

    Log("Building final MSA with %u sequences.\n", (unsigned)finalSeqs.size());

    // Create a new MSA from these labels & sequences
    MSA mCloaked;
    // We must strip the '>' from labels, because MSA::FromStrings2 expects them *without* a '>' prefix.
    for (auto& lab : finalLabels)
    {
        if (!lab.empty() && lab[0] == '>')
            lab.erase(0, 1);
    }
    mCloaked.FromStrings2(finalLabels, finalSeqs);

    mCloaked.ToFASTAFile(outputFilename);

    Log("Cloak: Output written to %s\n", outputFilename.c_str());
}