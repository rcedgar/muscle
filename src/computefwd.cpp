#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

/////////////////////////////////////////////////////////////////
// Forward probability matrix for pair of sequences.
// For efficiency, matrix is represented as an 
// array with the following indexing scheme:
//
//    forward[i + HMMSTATE_COUNT * (j * (L2+1) + k)]
//    = probability of aligning through j characters
//    of the first sequence, k characters of the second sequence,
//    and ending in state i.
/////////////////////////////////////////////////////////////////
vector<float> *PairHMM::ComputeForwardMatrix(Sequence* seq1, Sequence* seq2)
	{
	asserta(seq1 != 0);
	asserta(seq2 != 0);

	const int L1 = seq1->GetLength();
	const int L2 = seq2->GetLength();
	asserta(L1 > 0);
	asserta(L2 > 0);

	const byte *CharPtr1 = (const byte *) seq1->GetCharPtr1();
	const byte *CharPtr2 = (const byte *) seq2->GetCharPtr1();

	vector<float>* forwardPtr = new vector<float>(HMMSTATE_COUNT * (L1 + 1) * (L2 + 1), LOG_ZERO);
	assert(forwardPtr);
	vector<float>& forward = *forwardPtr;

	forward[0 + HMMSTATE_COUNT * (1 * (L2 + 1) + 1)] =
		m_StartScore[0] + m_MatchScore[CharPtr1[1]][CharPtr2[1]];

	for (int k = 0; k < InsertStateCount; k++)
		{
		forward[2 * k + 1 + HMMSTATE_COUNT * (1 * (L2 + 1) + 0)] =
			m_StartScore[2 * k + 1] + m_InsScore[CharPtr1[1]];

		forward[2 * k + 2 + HMMSTATE_COUNT * (0 * (L2 + 1) + 1)] =
			m_StartScore[2 * k + 2] + m_InsScore[CharPtr2[1]];
		}

// Offset for each index combination
	int ij = 0;
	int i1j = -L2 - 1;
	int ij1 = -1;
	int i1j1 = -L2 - 2;

	ij *= HMMSTATE_COUNT;
	i1j *= HMMSTATE_COUNT;
	ij1 *= HMMSTATE_COUNT;
	i1j1 *= HMMSTATE_COUNT;

// Forward scores
	for (int i = 0; i <= L1; i++)
		{
		byte c1 = CharPtr1[i];
		for (int j = 0; j <= L2; j++)
			{
			byte c2 = CharPtr2[j];

			if (i > 1 || j > 1)
				{
				if (i > 0 && j > 0)
					{
					forward[0 + ij] = forward[0 + i1j1] + m_TransScore[0][0];
					for (int k = 1; k < HMMSTATE_COUNT; k++)
						LOG_PLUS_EQUALS(forward[0 + ij], forward[k + i1j1] + m_TransScore[k][0]);
					forward[0 + ij] += m_MatchScore[c1][c2];
					}
				if (i > 0)
					{
					for (int k = 0; k < InsertStateCount; k++)
						forward[2 * k + 1 + ij] = m_InsScore[c1] +
						LOG_ADD(forward[0 + i1j] + m_TransScore[0][2 * k + 1],
							forward[2 * k + 1 + i1j] + m_TransScore[2 * k + 1][2 * k + 1]);
					}
				if (j > 0)
					{
					for (int k = 0; k < InsertStateCount; k++)
						forward[2 * k + 2 + ij] = m_InsScore[c2] +
						LOG_ADD(forward[0 + ij1] + m_TransScore[0][2 * k + 2],
							forward[2 * k + 2 + ij1] + m_TransScore[2 * k + 2][2 * k + 2]);
					}
				}

			ij += HMMSTATE_COUNT;
			i1j += HMMSTATE_COUNT;
			ij1 += HMMSTATE_COUNT;
			i1j1 += HMMSTATE_COUNT;
			}
		}

	return forwardPtr;
	}
