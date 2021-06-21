#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

/////////////////////////////////////////////////////////////////
// Backward probability matrix for pair of sequences.
// For efficiency, matrix is represented as an 
// array with the following indexing scheme:
//
//    backward[i + HMMSTATE_COUNT * (j * (L2+1) + k)]
//    = probability of starting in state i and aligning from 
//		character j+1 to the end of the first seq
//		and from character k+1 to the end of the second seq.
/////////////////////////////////////////////////////////////////
vector<float> *PairHMM::ComputeBackwardMatrix(Sequence* seq1, Sequence* seq2)
	{
	asserta(seq1 != 0);
	asserta(seq2 != 0);

	const int L1 = seq1->GetLength();
	const int L2 = seq2->GetLength();
	const byte *CharPtr1 = (const byte *) seq1->GetCharPtr1();
	const byte *CharPtr2 = (const byte *) seq2->GetCharPtr1();

	vector<float>* backwardPtr = new vector<float>(HMMSTATE_COUNT * (L1 + 1) * (L2 + 1), LOG_ZERO);
	assert(backwardPtr);
	vector<float>& backward = *backwardPtr;

// initialization condition
	for (int k = 0; k < HMMSTATE_COUNT; k++)
		backward[HMMSTATE_COUNT * ((L1 + 1) * (L2 + 1) - 1) + k] = m_StartScore[k];

// Offset for each index combination
	int ij = (L1 + 1) * (L2 + 1) - 1;
	int i1j = ij + L2 + 1;
	int ij1 = ij + 1;
	int i1j1 = ij + L2 + 2;

	ij *= HMMSTATE_COUNT;
	i1j *= HMMSTATE_COUNT;
	ij1 *= HMMSTATE_COUNT;
	i1j1 *= HMMSTATE_COUNT;

// Backward scores
	for (int i = L1; i >= 0; i--)
		{
		byte c1 = CharPtr1[i + 1];
		for (int j = L2; j >= 0; j--)
			{
			byte c2 = (j == L2) ? '~' : CharPtr2[j + 1];

			if (i < L1 && j < L2)
				{
				const float ProbXY = backward[0 + i1j1] + m_MatchScore[c1][c2];
				for (int k = 0; k < HMMSTATE_COUNT; k++)
					LOG_PLUS_EQUALS(backward[k + ij], ProbXY + m_TransScore[k][0]);
				}
			if (i < L1)
				{
				for (int k = 0; k < InsertStateCount; k++)
					{
					LOG_PLUS_EQUALS(backward[0 + ij], backward[2 * k + 1 + i1j] + m_InsScore[c1] + m_TransScore[0][2 * k + 1]);
					LOG_PLUS_EQUALS(backward[2 * k + 1 + ij], backward[2 * k + 1 + i1j] + m_InsScore[c1] + m_TransScore[2 * k + 1][2 * k + 1]);
					}
				}
			if (j < L2)
				{
				for (int k = 0; k < InsertStateCount; k++)
					{
					LOG_PLUS_EQUALS(backward[0 + ij], backward[2 * k + 2 + ij1] + m_InsScore[c2] + m_TransScore[0][2 * k + 2]);
					LOG_PLUS_EQUALS(backward[2 * k + 2 + ij], backward[2 * k + 2 + ij1] + m_InsScore[c2] + m_TransScore[2 * k + 2][2 * k + 2]);
					}
				}

			ij -= HMMSTATE_COUNT;
			i1j -= HMMSTATE_COUNT;
			ij1 -= HMMSTATE_COUNT;
			i1j1 -= HMMSTATE_COUNT;
			}
		}

	return backwardPtr;
	}
