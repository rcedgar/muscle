#include "muscle.h"
#include "best3.h"

float CalcAlnScoreFlat(const float *Post, uint LX, uint LY, float *DPRows)
	{
	float *Row = DPRows;
	for (uint j = 0; j <= LY; ++j)
		Row[j] = 0;

	const float *PostPtr = Post;
	for (uint i = 1; i <= LX; ++i)
		{
		float Currj1 = 0;
		float Prevj1 = Row[0];
		Row[0] = 0;
		for (uint j = 1; j <= LY; ++j)
			{
			float Prevj = Row[j];

			float P = *PostPtr++;
			float B = Prevj1 + P;
			float X = Prevj;
			float Y = Currj1;

			Prevj1 = Row[j];
			Best3(B, X, Y, &Currj1);
			Row[j] = Currj1;
			}
		}
	float Score = Row[LY];
	return Score;
	}
