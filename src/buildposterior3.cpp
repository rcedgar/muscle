#include "myutils.h"
#include "probcons.h"
#include "pairhmm.h"

#define TRACE	0

// Builds a posterior probability matrix needed to align a pair
// of alignments.  Mathematically, the returned matrix M is
// defined as follows:
//    M[SeqIndex1,SeqIndex2] =     sum          sum      f(s,t,SeqIndex1,SeqIndex2)
//                                  s in align1  t in align2
// where
//                                [  P(s[SeqIndex1] <--> t[SeqIndex2])
//                                [       if s[SeqIndex1] is a letter in the ith column of align1 and
//                                [          t[SeqIndex2] it a letter in the jth column of align2
//  f(s,t,SeqIndex1,SeqIndex2) =  [
//                                [  0    otherwise
//
// This is a variant of BuildPosterior() where sparseMatrices contains
// all pairs with one sequence from align1 and the other from align2,
// rather than all pairs in the union as in BuildPosterior.
vector<float>* PairHMM::BuildPosterior3(
  const MultiSequence* align1,
  const MultiSequence* align2,
  const vector<int> &SeqIndexes1,
  const vector<int> &SeqIndexes2,
  const vector<SparseMatrix*>& sparseMatrices)
	{
	const int SeqCount1 = align1->GetNumSequences();
	const int SeqCount2 = align2->GetNumSequences();

	const int ColCount1 = align1->GetSequence(0)->GetLength();
	const int ColCount2 = align2->GetSequence(0)->GetLength();

	const int PosteriorRowCount = (ColCount1 + 1);
	const int PosteriorColCount = (ColCount2 + 1);

	const int PosteriorSize = PosteriorRowCount*PosteriorColCount;
	const uint PairCount = SIZE(sparseMatrices);
	asserta(SIZE(SeqIndexes1) == PairCount);
	asserta(SIZE(SeqIndexes2) == PairCount);

	vector<float>* posteriorPtr = new vector<float>(PosteriorSize, 0);
	assert(posteriorPtr);
	vector<float>& posterior = *posteriorPtr;

	vector<uint> PosToCol1;
	vector<uint> PosToCol2;
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		uint SeqIndex1 = SeqIndexes1[PairIndex];
		uint SeqIndex2 = SeqIndexes2[PairIndex];

		const Sequence *Seq1 = align1->GetSequence(SeqIndex1);
		const Sequence *Seq2 = align2->GetSequence(SeqIndex2);

		const int ColCountSeq1 = Seq1->GetLength();
		const int ColCountSeq2 = Seq2->GetLength();

		asserta(ColCountSeq1 == ColCount1);
		asserta(ColCountSeq2 == ColCount2);

		Seq1->GetPosToCol_OneBased(PosToCol1);
		Seq2->GetPosToCol_OneBased(PosToCol2);

		const SparseMatrix* matrix = sparseMatrices[PairIndex];
		const int L1 = matrix->GetSeq1Length();
		const int L2 = matrix->GetSeq2Length();
		for (int Pos1 = 1; Pos1 <= L1; ++Pos1)
			{
			vector<PIF>::iterator row = matrix->GetRowPtr(Pos1);
			int Col1 = PosToCol1[Pos1];
			int base = Col1 * PosteriorColCount;
			int rowSize = matrix->GetRowSize(Pos1);

			for (int SparseCell2 = 0; SparseCell2 < rowSize; ++SparseCell2)
				{
				int Pos2 = row[SparseCell2].first;
				float Prob = row[SparseCell2].second;
				int Col2 = PosToCol2[Pos2];
				int Offset = base + Col2;
				asserta(Offset < PosteriorSize);
				posterior[Offset] += Prob;
				}
			}
		}
	return posteriorPtr;
	}
