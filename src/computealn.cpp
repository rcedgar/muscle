#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

// Computes alignment based on the given posterior matrix.
// This is done by finding the maximum summing path (or
// maximum weight trace) through the posterior matrix.
// The alignment is returned as a pair:
//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and Y's
//        denote insertions in one of the two sequences and B's
//        denote that both sequences are present (i.e. matches).
//    (2) a float indicating the sum achieved
// Traceback uses L, U, D (Left, Up, Diagonal) for X, Y, B.
pair<vector<char>*, float> PairHMM::ComputeAlignment(int seq1Length, int seq2Length,
  const vector<float>& posterior)
	{
	float* twoRows = new float[(seq2Length + 1) * 2];
	assert(twoRows);
	float* oldRow = twoRows;
	float* newRow = twoRows + seq2Length + 1;

	char* tracebackMatrix = new char[(seq1Length + 1) * (seq2Length + 1)];
	assert(tracebackMatrix);
	char* tracebackPtr = tracebackMatrix;

	vector<float>::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

	// initialization
	for (int i = 0; i <= seq2Length; i++)
		{
		oldRow[i] = 0;
		*(tracebackPtr++) = 'L';
		}

	// fill in matrix
	for (int i = 1; i <= seq1Length; i++)
		{
		// initialize left column
		newRow[0] = 0;
		posteriorPtr++;
		*(tracebackPtr++) = 'U';

		for (int j = 1; j <= seq2Length; j++)
			{
			ChooseBestOfThree(*(posteriorPtr++) + oldRow[j - 1],
				newRow[j - 1], oldRow[j],
				'D', 'L', 'U', &newRow[j], tracebackPtr++);
			}

		float* temp = oldRow;
		oldRow = newRow;
		newRow = temp;
		}

	float total = oldRow[seq2Length];
	delete[] twoRows;

	vector<char>* alignment = new vector<char>;
	assert(alignment);
	int r = seq1Length;
	int c = seq2Length;
	while (r != 0 || c != 0)
		{
		char ch = tracebackMatrix[r * (seq2Length + 1) + c];
		switch (ch)
			{
			case 'L': c--; alignment->push_back('Y'); break;
			case 'U': r--; alignment->push_back('X'); break;
			case 'D': c--; r--; alignment->push_back('B'); break;
			default: assert(false);
			}
		}

	delete[] tracebackMatrix;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, total);
	}
