#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

#define TRACE	0

// Builds a posterior probability matrix needed to align a pair
// of alignments.  Mathematically, the returned matrix M is
// defined as follows:
//    M[i,j] =     sum          sum      f(s,t,i,j)
//             s in align1  t in align2
// where
//                  [  P(s[i] <--> t[j])
//                  [       if s[i] is a letter in the ith column of align1 and
//                  [          t[j] it a letter in the jth column of align2
//    f(s,t,i,j) =  [
//                  [  0    otherwise
//
vector<float>* PairHMM::BuildPosterior5(MultiSequence* align1, MultiSequence* align2,
  const vector<vector<SparseMatrix*> > &sparseMatrices,
  const vector<vector<SparseMx> > &sparseMxs) const
	{
	const int seq1Length = align1->GetSequence(0)->GetLength();
	const int seq2Length = align2->GetSequence(0)->GetLength();

	vector<float>* posteriorPtr = 
	  new vector<float>((seq1Length + 1) * (seq2Length + 1), 0);
	assert(posteriorPtr);
	vector<float>& posterior = *posteriorPtr;

	// for each s in align1
	for (int i = 0; i < align1->GetNumSequences(); i++)
		{
		int first = align1->GetSequence(i)->GetIndex();
		vector<int>* mapping1 = align1->GetSequence(i)->GetPosToCol();

		// for each t in align2
		for (int j = 0; j < align2->GetNumSequences(); j++)
			{
			int second = align2->GetSequence(j)->GetIndex();
			vector<int>* mapping2 = align2->GetSequence(j)->GetPosToCol();

			if (first < second)
				{
				// get the associated sparse matrix
				SparseMatrix* matrix = sparseMatrices[first][second];

				for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++)
					{
					vector<PIF>::iterator row = matrix->GetRowPtr(ii);
					int base = (*mapping1)[ii] * (seq2Length + 1);
					int rowSize = matrix->GetRowSize(ii);

					// add in all relevant values
					for (int jj = 0; jj < rowSize; jj++)
						posterior[base + (*mapping2)[row[jj].first]] += row[jj].second;
					}
				}
			else
				{
				// get the associated sparse matrix
				SparseMatrix* matrix = sparseMatrices[second][first];

				for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++)
					{
					vector<PIF>::iterator row = matrix->GetRowPtr(jj);
					int base = (*mapping2)[jj];
					int rowSize = matrix->GetRowSize(jj);

					// add in all relevant values
					for (int ii = 0; ii < rowSize; ii++)
						posterior[base + (*mapping1)[row[ii].first] * (seq2Length + 1)] += row[ii].second;
					}
				}
			delete mapping2;
			}
		delete mapping1;
		}
#if TRACE
	{
	Log("BuildPosterior:\n");
	LogPosterior(*posteriorPtr, seq1Length, seq2Length);
	SparseMatrix *SM = new SparseMatrix(seq1Length, seq2Length, *posteriorPtr);
	SM->LogMe();
	delete SM;
	}
#endif
	return posteriorPtr;
	}
