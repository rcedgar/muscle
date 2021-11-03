#include "muscle.h"

float AlignPairFlat_SparsePost(const Sequence *Seq1, const Sequence *Seq2,
  string &Path, MySparseMx *SparsePost)
	{
	InitProbcons();

	uint L1 = Seq1->GetLength();
	uint L2 = Seq2->GetLength();
	asserta(L1 > 0);
	asserta(L2 > 0);

	const byte *ByteSeq1 = Seq1->GetBytePtr();
	const byte *ByteSeq2 = Seq2->GetBytePtr();

	float *Fwd = AllocFB(L1, L2);
	float *Bwd = AllocFB(L1, L2);
	float *Post = AllocPost(L1, L2);

	CalcFwdFlat(ByteSeq1, L1, ByteSeq2, L2, Fwd);
	CalcBwdFlat(ByteSeq1, L1, ByteSeq2, L2, Bwd);
	CalcPostFlat(Fwd, Bwd, L1, L2, Post);
	delete Fwd;
	delete Bwd;

	float *DPRows = AllocDPRows(L1, L2);
	char *TB = AllocTB(L1, L2);
	float Score = CalcAlnFlat(Post, L1, L2, DPRows, TB, Path);
	if (SparsePost != 0)
		SparsePost->FromPost(Post, L1, L2);
	delete Post;
	delete DPRows;
	delete TB;

	asserta(L1 > 0 && L2 > 0);
	float EA = Score/min(L1, L2);
	return EA;
	}

float AlignPairFlat(const Sequence *Seq1, const Sequence *Seq2, string &Path)
	{
	float EA = AlignPairFlat_SparsePost(Seq1, Seq2, Path, 0);
	return EA;
	}
