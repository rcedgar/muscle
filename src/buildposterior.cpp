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
//                  [          t[j] is a letter in the jth column of align2
//    f(s,t,i,j) =  [
//                  [  0    otherwise
//
vector<float>* PairHMM::BuildPosterior(MultiSequence* align1,
  MultiSequence* align2,
  const vector<vector<SparseMatrix*> >& sparseMatrices)
	{
	const int seq1Length = align1->GetSequence(0)->GetLength();
	const int seq2Length = align2->GetSequence(0)->GetLength();

	const uint SMSize = SIZE(sparseMatrices);
	asserta(SMSize > 0);
	asserta(SIZE(sparseMatrices[0]) == SMSize);

	vector<float>* posteriorPtr = 
	  new vector<float>((seq1Length + 1) * (seq2Length + 1), 0);
	assert(posteriorPtr);
	vector<float>& posterior = *posteriorPtr;

// for each s in align1
	vector<uint> PosToCol1;
	vector<uint> PosToCol2;
	for (int i = 0; i < align1->GetNumSequences(); i++)
		{
		uint SMI_i = align1->GetSequence(i)->GetSMI();
		asserta(SMI_i != UINT_MAX);
		asserta(SMI_i < SMSize);

		align1->GetSequence(i)->GetPosToCol_OneBased(PosToCol1);

	// for each t in align2
		for (int j = 0; j < align2->GetNumSequences(); j++)
			{
			uint SMI_j = align2->GetSequence(j)->GetSMI();
			asserta(SMI_j != UINT_MAX);
			asserta(SMI_j < SMSize);

			align2->GetSequence(j)->GetPosToCol_OneBased(PosToCol2);

			if (SMI_i < SMI_j)
				{
				SparseMatrix* matrix = sparseMatrices[SMI_i][SMI_j];
				asserta(matrix != 0);

				for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++)
					{
					vector<PIF>::iterator row = matrix->GetRowPtr(ii);
					int base = PosToCol1[ii] * (seq2Length + 1);
					int rowSize = matrix->GetRowSize(ii);

					for (int jj = 0; jj < rowSize; jj++)
						posterior[base + PosToCol2[row[jj].first]] += row[jj].second;
					}
				}
			else
				{
				SparseMatrix* matrix = sparseMatrices[SMI_j][SMI_i];
				asserta(matrix != 0);

				for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++)
					{
					vector<PIF>::iterator row = matrix->GetRowPtr(jj);
					int base = PosToCol2[jj];
					int rowSize = matrix->GetRowSize(jj);

					for (int ii = 0; ii < rowSize; ii++)
						posterior[base + PosToCol1[row[ii].first] * (seq2Length + 1)] += row[ii].second;
					}
				}
			}
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
