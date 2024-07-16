#include "muscle.h"
#include "mega.h"

/***
Fwd[s][i][j] = 
	probability of aligning 
	first i letters of X to
	first j letters of Y and
	ending in state s.
***/

void CalcFwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat)
	{
	asserta(!Mega::m_Loaded);

#include "hmmscores.h"
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);

	char x0 = X[0];
	char y0 = Y[0];
	float Ins_x0 = InsScore[x0];
	float Ins_y0 = InsScore[y0];
	float Emit_x0_y0 = MatchScore[x0][y0];

	const uint LY1 = LY+1;
	const uint Base_0_0 = HMMSTATE_COUNT*(0*(LY1) + 0);
	const uint Base_1_1 = HMMSTATE_COUNT*(1*(LY1) + 1);
	const uint Base_1_0 = HMMSTATE_COUNT*(1*(LY1) + 0);
	const uint Base_0_1 = HMMSTATE_COUNT*(0*(LY1) + 1);

	const uint BaseInc_i = HMMSTATE_COUNT*LY1;
	const uint BaseInc_j = HMMSTATE_COUNT;

	Flat[Base_0_0 + HMMSTATE_M] = LOG_ZERO;		// M(0,0)
	Flat[Base_0_0 + HMMSTATE_IX] = LOG_ZERO;	// IX(0,0)
	Flat[Base_0_0 + HMMSTATE_JX] = LOG_ZERO;	// JX(0,0)
	Flat[Base_0_0 + HMMSTATE_IY] = LOG_ZERO;	// IY(0,0)
	Flat[Base_0_0 + HMMSTATE_JY] = LOG_ZERO;	// JY(0,0)

	Flat[Base_1_1 + HMMSTATE_M] = tSM + Emit_x0_y0;
	Flat[Base_1_0 + HMMSTATE_IX] = tSI + Ins_x0;
	Flat[Base_1_0 + HMMSTATE_JX] = tSJ + Ins_x0;
	Flat[Base_0_1 + HMMSTATE_IY] = tSI + Ins_y0;
	Flat[Base_0_1 + HMMSTATE_JY] = tSJ + Ins_y0;

	uint Base = Base_1_0;
	for (uint i = 1; i <= LX; ++i)
		{
		Flat[Base + HMMSTATE_M] = LOG_ZERO;
		Flat[Base + HMMSTATE_IY] = LOG_ZERO;
		Flat[Base + HMMSTATE_JY] = LOG_ZERO;

		Base += BaseInc_i;
		}

	Base = Base_0_1;
	for (uint j = 1; j <= LY; ++j)
		{
		Flat[Base + HMMSTATE_M] = LOG_ZERO;
		Flat[Base + HMMSTATE_IX] = LOG_ZERO;
		Flat[Base + HMMSTATE_JX] = LOG_ZERO;

		Base += BaseInc_j;
		}

	Base = Base_1_0;
	uint NextBase = Base + BaseInc_i;
	for (uint i = 1; i < LX; ++i)
		{
		char x = X[i];
		float Emit_x = InsScore[x];

		Flat[NextBase + HMMSTATE_IX] = Flat[Base + HMMSTATE_IX] + tII + Emit_x;
		Flat[NextBase + HMMSTATE_JX] = Flat[Base + HMMSTATE_JX] + tJJ + Emit_x;

		Base = NextBase;
		NextBase += BaseInc_i;
		}

	Base = Base_0_1;
	NextBase = Base + BaseInc_j;
	for (uint j = 1; j < LY; ++j)
		{
		char y = Y[j];
		float Emit_y = InsScore[y];

		Flat[NextBase + HMMSTATE_IY] = Flat[Base + HMMSTATE_IY] + tII + Emit_y;
		Flat[NextBase + HMMSTATE_JY] = Flat[Base + HMMSTATE_JY] + tJJ + Emit_y;

		Base = NextBase;
		NextBase += BaseInc_j;
		}

	uint Base_i_j = Base_1_1;
	uint Base_i1_j = Base_0_1;
	uint Base_i_j1 = Base_1_0;
	uint Base_i1_j1 = Base_0_0;

	for (uint i = 1; i <= LX; ++i)
		{
		char x = X[i-1];
		float Emit_x = InsScore[x];

		for (uint j = 1; j <= LY; ++j)
			{
			char y = Y[j-1];
			float Emit_y = InsScore[y];
			float Emit_Pair = MatchScore[x][y];
			if (i == 1 && j == 1)
				Flat[Base_1_1 + HMMSTATE_M] = tSM + Emit_x0_y0;
			else
				{
				float M_M = Flat[Base_i1_j1 + HMMSTATE_M] + tMM;
				float IX_M = Flat[Base_i1_j1 + HMMSTATE_IX] + tIM;
				float JX_M = Flat[Base_i1_j1 + HMMSTATE_JX] + tJM;
				float IY_M = Flat[Base_i1_j1 + HMMSTATE_IY] + tIM;
				float JY_M = Flat[Base_i1_j1 + HMMSTATE_JY] + tJM;

				float SumPrev = LOG_ADD(M_M, IX_M, JX_M, IY_M, JY_M);
				Flat[Base_i_j + HMMSTATE_M] = SumPrev + Emit_Pair;
				}

			float PrevM_i1_j = Flat[Base_i1_j + HMMSTATE_M];
			float PrevM_i_j1 = Flat[Base_i_j1 + HMMSTATE_M];
			float M_IX = PrevM_i1_j + tMI;
			float IX_IX = Flat[Base_i1_j + HMMSTATE_IX] + tII;
			Flat[Base_i_j + HMMSTATE_IX] = LOG_ADD(IX_IX,  M_IX) + Emit_x;

			float M_JX = PrevM_i1_j + tMJ;
			float JX_JX = Flat[Base_i1_j + HMMSTATE_JX] + tJJ;
			Flat[Base_i_j + HMMSTATE_JX] = LOG_ADD(JX_JX, M_JX) + Emit_x;

			float M_IY = PrevM_i_j1 + tMI;
			float IY_IY = Flat[Base_i_j1 + HMMSTATE_IY] + tII;
			Flat[Base_i_j + HMMSTATE_IY] = LOG_ADD(IY_IY, M_IY) + Emit_y;

			float M_JY = PrevM_i_j1 + tMJ;
			float JY_JY = Flat[Base_i_j1 + HMMSTATE_JY] + tJJ;
			Flat[Base_i_j + HMMSTATE_JY] = LOG_ADD(JY_JY, M_JY) + Emit_y;

			Base_i_j += BaseInc_j;
			Base_i1_j += BaseInc_j;
			Base_i_j1 += BaseInc_j;
			Base_i1_j1 += BaseInc_j;
			}

		Base_i_j += BaseInc_j;
		Base_i1_j += BaseInc_j;
		Base_i_j1 += BaseInc_j;
		Base_i1_j1 += BaseInc_j;
		}
	}
