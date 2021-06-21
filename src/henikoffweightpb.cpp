#include "muscle.h"
#include "msa.h"

/***
Compute Henikoff weights.
Steven Henikoff and Jorja G. Henikoff (1994), Position-based sequence weights.
J. Mol. Biol., 243(4):574-578.

Award each different residue an equal share of the weight, and then to divide up
that weight equally among the sequences sharing the same residue. So if in a
position of a multiple alignment, r different residues are represented, a residue
represented in only one sequence contributes a score of 1/r to that sequence, whereas a
residue represented in s sequences contributes a score of 1/rs to each of the s
sequences. For each sequence, the contributions from each position are summed to give
a sequence weight.

Here we use the variant from PSI-BLAST, which (a) treats gaps as a 21st letter,
and (b) ignores columns that are perfectly conserved.

>>> WARNING -- I SUSPECT THIS DOESN'T WORK CORRECTLY <<<
***/

void MSA::CalcHenikoffWeightsColPB(unsigned uColIndex) const
	{
	const unsigned uSeqCount = GetSeqCount();

// Compute letter counts in this column
	unsigned uLetterCount[MAX_ALPHA+1];
	memset(uLetterCount, 0, (MAX_ALPHA+1)*sizeof(unsigned));
	unsigned uLetter;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		if (IsGap(uSeqIndex, uColIndex) || IsWildcard(uSeqIndex, uColIndex))
			uLetter = MAX_ALPHA;
		else
			uLetter = GetLetter(uSeqIndex, uColIndex);
		++(uLetterCount[uLetter]);
		}

// Check for special case of perfect conservation
	for (unsigned uLetter = 0; uLetter < MAX_ALPHA+1; ++uLetter)
		{
		unsigned uCount = uLetterCount[uLetter];
		if (uCount > 0)
			{
		// Perfectly conserved?
			if (uCount == uSeqCount)
				return;
			else
			// If count > 0 but less than nr. sequences, can't be conserved
				break;
			}
		}

// Compute weight contributions
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uLetter;
		if (IsGap(uSeqIndex, uColIndex) || IsWildcard(uSeqIndex, uColIndex))
			uLetter = MAX_ALPHA;
		else
			uLetter = GetLetter(uSeqIndex, uColIndex);
		const unsigned uCount = uLetterCount[uLetter];
		m_Weights[uSeqIndex] += (WEIGHT) (1.0/uCount);
		}
	}

bool MSA::IsGapSeq(unsigned uSeqIndex) const
	{
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		if (!IsGap(uSeqIndex, uColIndex))
			return false;
	return true;
	}

void MSA::SetUniformWeights() const
	{
	const unsigned uSeqCount = GetSeqCount();
	if (0 == uSeqCount)
		return;

	const WEIGHT w = (WEIGHT) (1.0 / uSeqCount);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		m_Weights[uSeqIndex] = w;
	}

void MSA::SetHenikoffWeightsPB() const
	{
	const unsigned uColCount = GetColCount();
	const unsigned uSeqCount = GetSeqCount();

	if (0 == uSeqCount)
		return;
	else if (1 == uSeqCount)
		{
		m_Weights[0] = 1.0;
		return;
		}
	else if (2 == uSeqCount)
		{
		m_Weights[0] = (WEIGHT) 0.5;
		m_Weights[1] = (WEIGHT) 0.5;
		return;
		}

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		m_Weights[uSeqIndex] = 0.0;

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		CalcHenikoffWeightsColPB(uColIndex);

// Set all-gap seqs weight to 0
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		if (IsGapSeq(uSeqIndex))
			m_Weights[uSeqIndex] = 0.0;

// Check for special case of identical sequences, which will cause all
// columns to be skipped becasue they're perfectly conserved.
	if (VectorIsZero(m_Weights, uSeqCount))
		VectorSet(m_Weights, uSeqCount, 1.0);

	Normalize(m_Weights, uSeqCount);
	}
