#include "muscle.h"

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
// This is a variant of BuildPosterior() where sparse posterior matrices 
// contain all pairs with one sequence from MSA1 and the other from MSA2,
// rather than all pairs in the union as in CalcPostFlat.

void CalcPosteriorFlat3(const MultiSequence &MSA1,
  const MultiSequence &MSA2,
  const vector<uint> &SeqIndexes1,
  const vector<uint> &SeqIndexes2,
  const vector<MySparseMx *> &SparseMxs,
  float *Flat)
	{
	const uint SeqCount1 = MSA1.GetSeqCount();
	const uint SeqCount2 = MSA1.GetSeqCount();

	const uint ColCount1 = MSA1.GetColCount();
	const uint ColCount2 = MSA2.GetColCount();

	const uint FlatSize = ColCount1*ColCount2;
	for (uint i = 0; i < FlatSize; ++i)
		Flat[i] = 0;

	vector<uint> PosToCol1;
	vector<uint> PosToCol2;

// May be subset of all pairs due to sampling
	const uint PairCount = SIZE(SparseMxs);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		uint SeqIndex1 = SeqIndexes1[PairIndex];
		uint SeqIndex2 = SeqIndexes2[PairIndex];

		const Sequence *Seq1 = MSA1.GetSequence(SeqIndex1);
		const Sequence *Seq2 = MSA2.GetSequence(SeqIndex2);

		const uint ColCountSeq1 = Seq1->GetLength();
		const uint ColCountSeq2 = Seq2->GetLength();

		asserta(ColCountSeq1 == ColCount1);
		asserta(ColCountSeq2 == ColCount2);

		const MySparseMx &PostMx12 = *SparseMxs[PairIndex];
		const uint L1 = PostMx12.GetLX();
		const uint L2 = PostMx12.GetLY();

		Seq1->GetPosToCol(PosToCol1);
		Seq2->GetPosToCol(PosToCol2);

		asserta(SIZE(PosToCol1) == L1);
		asserta(SIZE(PosToCol2) == L2);

		for (uint Pos1 = 0; Pos1 < L1; ++Pos1)
			{
			uint Offset = PostMx12.GetOffset(Pos1);
			uint RowSize = PostMx12.GetSize(Pos1);
			assert(Pos1 < SIZE(PosToCol1));
			uint Col1 = PosToCol1[Pos1];
			uint FlatBase = Col1 * ColCount2;

			for (uint k = 0; k < RowSize; ++k)
				{
				float Prob = PostMx12.GetProb_Offset(Offset + k);
				uint Pos2 = PostMx12.GetCol_Offset(Offset + k);
				assert(Pos2 < SIZE(PosToCol2));
				uint Col2 = PosToCol2[Pos2];
				uint FlatOffset = FlatBase + Col2;
				assert(FlatOffset < FlatSize);
				Flat[FlatOffset] += Prob;
				}
			}
		}
	}
