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
 *    - Columns not meeting the minimum occupancy threshold will be discarded. -mincol 2, the
 *      default, is the published algorithm; higher values are an extension.
 *    - Every input sequence is written, in label order, even if it keeps no residue.
 *    - Gaps are '-' and '.', as in MUSCLE's loader; output gaps are '-'.
 *    - Every error, including no surviving column, a failed write and running out of memory,
 *      exits 1 with a message on stderr and leaves no output file.
 *    - -output - writes the result to standard output.
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
#include <cerrno>
#include <cstring>
#include <new>
#include <sys/stat.h>


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
//   (starting at 1) to each residue; every gap ('-' or '.', MUSCLE's isgap) becomes '-'.
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
                if (isgap(letter)) {
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


// Builds the numeric CLOAK output: one entry per retained column, each a vector over rows that
// holds a residue number (1-based, as made by convert_to_nums) or '-'.
// Two residues of a column of alignment 0 are kept together iff they share a column in every
// alignment. Each residue is keyed by the columns it occupies in alignments 1..K-1, so the
// retained columns are exactly the groups of residues with equal keys. Groups are emitted in
// lexicographic key order, which is the order in which the recursive search used up to f3559af
// first reached them, so the output file is unchanged; that search reached each group many times
// over and stored every copy, which exhausted memory on ensembles whose alignments disagree.
// Groups with fewer than minimumColumns residues are dropped here, so remove_singletons, which
// filters again after conversion to letters, removes nothing.
std::vector<std::vector<NumOrDash>> create_output(
    const std::vector<std::vector<std::vector<NumOrDash>>>& all_als,
    uint minimumColumns)
{
    std::vector<std::vector<NumOrDash>> output;
    const size_t K = all_als.size();
    const auto& first_alignment = all_als[0];
    const size_t n = first_alignment.size();
    if (n == 0 || first_alignment[0].empty())
        return output;

    // pos_col[k][r][p-1] = column of residue number p of row r in alignment k (k >= 1).
    std::vector<size_t> residues0(n, 0);
    for (size_t r = 0; r < n; ++r)
        for (const auto& v : first_alignment[r])
            if (v.type == NumOrDashType::INTEGER)
                ++residues0[r];
    std::vector<std::vector<std::vector<int>>> pos_col(K);
    for (size_t k = 1; k < K; ++k) {
        pos_col[k].resize(n);
        for (size_t r = 0; r < n; ++r) {
            const auto& row = all_als[k][r];
            for (size_t c = 0; c < row.size(); ++c)
                if (row[c].type == NumOrDashType::INTEGER)
                    pos_col[k][r].push_back(static_cast<int>(c));
            // The loader accepts only MSAs whose ungapped sequences agree, so counts are equal.
            asserta(pos_col[k][r].size() == residues0[r]);
        }
    }

    const size_t num_cols = first_alignment[0].size();
    std::vector<std::pair<std::vector<int>, size_t>> keyed;   // (key, row)
    for (size_t x = 0; x < num_cols; ++x) {
        keyed.clear();
        for (size_t r = 0; r < n; ++r) {
            const NumOrDash& v = first_alignment[r][x];
            if (v.type != NumOrDashType::INTEGER)
                continue;
            const size_t p = static_cast<size_t>(v.intValue - 1);
            std::vector<int> key(K - 1);
            for (size_t k = 1; k < K; ++k)
                key[k - 1] = pos_col[k][r][p];
            keyed.emplace_back(std::move(key), r);
        }
        std::sort(keyed.begin(), keyed.end());
        for (size_t i = 0; i < keyed.size(); ) {
            size_t j = i + 1;
            while (j < keyed.size() && keyed[j].first == keyed[i].first)
                ++j;
            if (j - i >= minimumColumns) {
                std::vector<NumOrDash> column(n, NumOrDash('-'));
                for (size_t t = i; t < j; ++t)
                    column[keyed[t].second] = first_alignment[keyed[t].second][x];
                output.push_back(std::move(column));
            }
            i = j;
        }
    }
    return output;
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


// Takes the output columns (each a vector over rows, gaps '-') and returns only those with at
// least minColumns residues.
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

    // Build an "intermediary" 2D array that collects the residues (non-gap
    // characters) column by column from first_alignment.
    size_t num_cols = first_alignment[0].size();
    for (size_t i = 0; i < num_cols; i++) {
        std::string s;
        // Concatenate the i-th element from every row
        for (const auto& row : first_alignment) {
            if (i < row.size()) {
                s.push_back(row[i]);
            }
        }
        // Build an array of characters from s, skipping gaps.
        std::vector<char> some_array;
        for (char ch : s) {
            if (isgap(ch)) {
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


// Writes the FASTA bytes of MSA::ToFASTAFile (">label", then lines of 60 characters, LF) with every
// write checked: MUSCLE's TextFile checks writes only with assert, so a full disk exited 0 with a
// truncated file. "-" is standard output, as for TextFile. A partial regular file is removed before
// Die, so a failed write leaves no output file.
static void write_fasta(const std::string& path, const std::vector<std::string>& labels,
    const std::vector<std::string>& seqs)
{
    const bool toStdout = (path == "-");
    FILE* f = toStdout ? stdout : fopen(path.c_str(), "wb");
    if (f == 0)
        Die("cloak: cannot create %s, errno=%d %s", path.c_str(), errno, strerror(errno));
    bool ok = true;
    for (size_t r = 0; ok && r < labels.size(); ++r) {
        ok = fprintf(f, ">%s\n", labels[r].c_str()) >= 0;
        const std::string& s = seqs[r];
        for (size_t i = 0; ok && i < s.size(); i += 60) {
            const size_t len = std::min<size_t>(60, s.size() - i);
            ok = fwrite(s.data() + i, 1, len, f) == len && fputc('\n', f) != EOF;
        }
    }
    ok = (toStdout ? fflush(f) == 0 : fclose(f) == 0) && ok;
    if (!ok) {
        struct stat st;
        if (!toStdout && stat(path.c_str(), &st) == 0 && (st.st_mode & S_IFMT) == S_IFREG)
            remove(path.c_str());
        Die("cloak: cannot write %s, errno=%d %s", path.c_str(), errno, strerror(errno));
    }
}


static void run_cloak()

{
    // Parse user options
    std::string ensemblePath = opt(cloak);  // user runs: muscle -cloak myfile ...
    if (ensemblePath.empty())
        Die("cloak: you must specify -cloak <ensembleFile or list-of-FASTA>");

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
        Die("cloak: -mincol must be >= 2");

    Log("Starting cloak. Ensemble path: %s, Output filename: %s, Minimum columns requested: %u\n",
        ensemblePath.c_str(), outputFilename.c_str(), userMincol);

    // Load the entire ensemble from the user-supplied file
    Ensemble E;
    E.FromFile(ensemblePath);
    if (E.GetMSACount() == 0)   // a path list naming no MSA loads with only a warning
        Die("cloak: the ensemble %s holds no MSA", ensemblePath.c_str());

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

    // 2) Create numeric output (each retained column once, so no duplicate removal is needed)
    Log("Creating initial numeric output.\n");
    auto nums_array_output = create_output(nums_answer, minimumColumns);
    if (nums_array_output.empty())
        Die("cloak: no column has %u or more residues", minimumColumns);

    // 3) Use the first alignment in three_d_array_letters as reference.
    auto reference_alignment = three_d_array_letters[0];

    // 4) Convert numeric output back to letters
    Log("Converting numeric output back to letters.\n");
    auto transposed_reference = transpose_output(reference_alignment);
    auto letters_array_output_initial =
        produce_output_in_letters(nums_array_output, transposed_reference);

    Log("Removing singletons with minimum column threshold = %u.\n", minimumColumns);
    auto letters_array_output = remove_singletons(letters_array_output_initial, minimumColumns);

    letters_array_output = transpose_output(letters_array_output);

    // 5) Every input sequence is written, under its label exactly as read and in label order,
    // however many residues it keeps (decision 007). allLabels[0] is in that order.
    const auto& firstLabels = allLabels[0];
    std::vector<std::string> finalSeqs;
    for (const auto& row : letters_array_output)
        finalSeqs.emplace_back(row.begin(), row.end());

    Log("Writing %u sequences.\n", (unsigned)finalSeqs.size());
    write_fasta(outputFilename, firstLabels, finalSeqs);

    Log("Cloak: Output written to %s\n", outputFilename.c_str());
}


void cmd_cloak()
{
    // Allocation failures exit 1 like every other error; uncaught, std::bad_alloc aborted (exit 134).
    // Unwinding frees run_cloak's memory before Die runs.
    try
    {
        run_cloak();
    }
    catch (const std::bad_alloc&)
    {
        Die("cloak: out of memory");
    }
    catch (...)
    {
        Die("cloak: unexpected exception");
    }
}