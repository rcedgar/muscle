#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

/////////////////////////////////////////////////////////////////
// PairHMM::ComputePosteriorMatrix()
//
// Computes the posterior probability matrix based on
// the forward and backward matrices.
/////////////////////////////////////////////////////////////////

vector<float>* PairHMM::ComputePosteriorMatrix(Sequence* seq1, Sequence* seq2,
  const vector<float>& forward, const vector<float>& backward)
	{
	assert(seq1);
	assert(seq2);

	const int seq1Length = seq1->GetLength();
	const int seq2Length = seq2->GetLength();

	float totalProb = ComputeTotalProbability(seq1Length, seq2Length,
		forward, backward);

	// compute posterior matrices
	vector<float>* posteriorPtr = new vector<float>((seq1Length + 1) * (seq2Length + 1)); assert(posteriorPtr);
	vector<float>& posterior = *posteriorPtr;

	int ij = 0;
	vector<float>::iterator ptr = posterior.begin();

	for (int i = 0; i <= seq1Length; i++)
		{
		for (int j = 0; j <= seq2Length; j++)
			{
			*(ptr++) = EXP(min(LOG_ONE, forward[ij] + backward[ij] - totalProb));
			ij += HMMSTATE_COUNT;
			}
		}

	posterior[0] = 0;

	return posteriorPtr;
	}
