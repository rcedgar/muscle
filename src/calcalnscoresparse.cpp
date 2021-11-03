#include "muscle.h"
#include "best3.h"

float CalcAlnScoreSparse(const MySparseMx &Mx)
	{
	uint LX = Mx.m_LX;
	uint LY = Mx.m_LY;

#if 0
	float Sum = 0;
	for (uint i = 0; i < Mx.m_LX; ++i)
		Sum += Mx.GetMaxProbRow(i);
	return Sum;
#endif

	float *Post = AllocPost(LX, LY);
	Mx.ToPost(Post);
	float *DPRows = AllocDPRows(LX, LY);
	float Score = CalcAlnScoreFlat(Post, LX, LY, DPRows);
	myfree(Post);
	myfree(DPRows);
	return Score;
	}
