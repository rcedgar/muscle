#include "muscle.h"

float AlignPairFlat_SparsePost(const string &Label1, const string &Label2,
  string &Path, MySparseMx *SparsePost)
	{
	uint L1 = GetSeqLengthByGlobalLabel(Label1);
	uint L2 = GetSeqLengthByGlobalLabel(Label2);
	float *Post = CalcPost(Label1, Label2);
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

float AlignPairFlat(const string &Label1, const string &Label2, string &Path)
	{
	float EA = AlignPairFlat_SparsePost(Label1, Label2, Path, 0);
	return EA;
	}
