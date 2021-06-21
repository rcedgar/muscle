#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

/////////////////////////////////////////////////////////////////
// PairHMM::ComputeTotalProbability()
//
// Computes the total probability of an alignment given
// the forward and backward matrices.
/////////////////////////////////////////////////////////////////

float PairHMM::ComputeTotalProbability(int seq1Length, int seq2Length,
  const vector<float>& forward, const vector<float>& backward)
	{
	// compute total probability
	float totalForwardProb = LOG_ZERO;
	float totalBackwardProb = LOG_ZERO;
	for (int k = 0; k < HMMSTATE_COUNT; k++) {
		LOG_PLUS_EQUALS(totalForwardProb,
			forward[k + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)] +
			backward[k + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)]);
		}

	totalBackwardProb =
		forward[0 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 1)] +
		backward[0 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 1)];

	for (int k = 0; k < InsertStateCount; k++) {
		LOG_PLUS_EQUALS(totalBackwardProb,
			forward[2 * k + 1 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 0)] +
			backward[2 * k + 1 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 0)]);
		LOG_PLUS_EQUALS(totalBackwardProb,
			forward[2 * k + 2 + HMMSTATE_COUNT * (0 * (seq2Length + 1) + 1)] +
			backward[2 * k + 2 + HMMSTATE_COUNT * (0 * (seq2Length + 1) + 1)]);
		}

	return (totalForwardProb + totalBackwardProb) / 2;
	}
