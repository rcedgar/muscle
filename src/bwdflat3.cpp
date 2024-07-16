#include "muscle.h"

/***
Bwd[s][i][j] = 
	probability of starting in state s and aligning
	last (LX-i) letters of X to
	last (LY-j) letters of Y.
***/

void CalcBwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat)
	{
	asserta(!Mega::m_Loaded);

#include "hmmscores.h"
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);

	const int iLX = int(LX);
	const int iLY = int(LY);

	const int LY1 = LY+1;
	const int BaseInc_i = HMMSTATE_COUNT*LY1;
	const int BaseInc_j = HMMSTATE_COUNT;

	uint Base = HMMSTATE_COUNT*(0*(LY1) + LY);
	for (int i = 0; i < iLX; ++i)
		{
		Flat[Base + HMMSTATE_IY] = LOG_ZERO;
		Flat[Base + HMMSTATE_JY] = LOG_ZERO;
		Base += BaseInc_i;
		}

	Base = HMMSTATE_COUNT*(LX*(LY1) + 0);
	for (int j = 0; j < iLY; ++j)
		{
		Flat[Base + HMMSTATE_IX] = LOG_ZERO;
		Flat[Base + HMMSTATE_JX] = LOG_ZERO;
		Base += BaseInc_j;
		}

	int Base_i_j = (int) HMMSTATE_COUNT*(LX*(LY1) + LY);
	int Base_i1_j = Base_i_j + BaseInc_i;
	int Base_i_j1 = Base_i_j + BaseInc_j;
	int Base_i1_j1 = Base_i_j + BaseInc_i + BaseInc_j;

	for (int i = iLX; i >= 0; --i)
		{
		char x = (i == iLX ? 0 : X[i]);
		float Emit_x = InsScore[x];

		for (int j = iLY; j >= 0; --j)
			{
			if (i == LX && j == LY)
				{
			// Special case for end-of-alignment
				Flat[Base_i_j + HMMSTATE_M] = tSM;
				Flat[Base_i_j + HMMSTATE_IX] = tSI;
				Flat[Base_i_j + HMMSTATE_IY] = tSI;
				Flat[Base_i_j + HMMSTATE_JX] = tSJ;
				Flat[Base_i_j + HMMSTATE_JY] = tSJ;

				Base_i_j -= BaseInc_j;
				Base_i1_j -= BaseInc_j;
				Base_i_j1 -= BaseInc_j;
				Base_i1_j1 -= BaseInc_j;
				continue;
				}

			char y = (j == iLY ? 0 : Y[j]);
			float Emit_y = InsScore[y];
			float Emit_xy = MatchScore[x][y];

			if (i < iLX && j < iLY)
				{
				float NextM  = Flat[Base_i1_j1 + HMMSTATE_M] + Emit_xy;
				float NextIX = Flat[Base_i1_j + HMMSTATE_IX] + Emit_x;
				float NextJX = Flat[Base_i1_j + HMMSTATE_JX] + Emit_x;
				float NextIY = Flat[Base_i_j1 + HMMSTATE_IY] + Emit_y;
				float NextJY = Flat[Base_i_j1 + HMMSTATE_JY] + Emit_y;

				if (i > 0 && j > 0)
					{
					float M_M  = tMM + NextM;
					float M_IX = tMI + NextIX;
					float M_JX = tMJ + NextJX;
					float M_IY = tMI + NextIY;
					float M_JY = tMJ + NextJY;
					Flat[Base_i_j + HMMSTATE_M] = LOG_ADD(M_M, M_IX, M_JX, M_IY, M_JY);
					}
				else
					Flat[Base_i_j + HMMSTATE_M] = LOG_ZERO;

				if (i > 0)
					{
					float IX_IX = tII + NextIX;
					float IX_M = tIM + NextM;
					Flat[Base_i_j + HMMSTATE_IX] = LOG_ADD(IX_IX,  IX_M);

					float JX_JX = tJJ + NextJX;
					float JX_M = tJM + NextM;
					Flat[Base_i_j + HMMSTATE_JX] = LOG_ADD(JX_JX,  JX_M);
					}
				else
					{
					Flat[Base_i_j + HMMSTATE_IX] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_JX] = LOG_ZERO;
					}

				if (j > 0)
					{
					float IY_IY = tII + NextIY;
					float IY_M = tIM + NextM;
					Flat[Base_i_j + HMMSTATE_IY] = LOG_ADD(IY_IY,  IY_M);

					float JY_JY = tJJ + NextJY;
					float JY_M = tJM + NextM;
					Flat[Base_i_j + HMMSTATE_JY] = LOG_ADD(JY_JY,  JY_M);
					}
				else
					{
					Flat[Base_i_j + HMMSTATE_IY] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_JY] = LOG_ZERO;
					}

				Base_i_j -= BaseInc_j;
				Base_i1_j -= BaseInc_j;
				Base_i_j1 -= BaseInc_j;
				Base_i1_j1 -= BaseInc_j;
				continue;
				}

			if (i < iLX)
				{
				assert(j == iLY);
				if (i > 0)
					{
					float NextIX = Flat[Base_i1_j + HMMSTATE_IX] + Emit_x;
					float NextJX = Flat[Base_i1_j + HMMSTATE_JX] + Emit_x;

					float M_IX = tMI + NextIX;
					float M_JX = tMJ + NextJX;

					Flat[Base_i_j + HMMSTATE_M] = LOG_ADD(M_IX, M_JX);
					Flat[Base_i_j + HMMSTATE_IX] = tII + NextIX;
					Flat[Base_i_j + HMMSTATE_JX] = tJJ + NextJX;
					}
				else
					{
					Flat[Base_i_j + HMMSTATE_M] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_IX] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_JX] = LOG_ZERO;
					}
				}

			if (j < iLY)
				{
				assert(i == iLX);

				float NextIY = Flat[Base_i_j1 + HMMSTATE_IY] + Emit_y;
				float NextJY = Flat[Base_i_j1 + HMMSTATE_JY] + Emit_y;

				float M_IY = tMI + NextIY;
				float M_JY = tMJ + NextJY;
				if (j > 0)
					{
					Flat[Base_i_j + HMMSTATE_M] = LOG_ADD(M_IY, M_JY);
					Flat[Base_i_j + HMMSTATE_IY] = tII + NextIY;
					Flat[Base_i_j + HMMSTATE_JY] = tJJ + NextJY;
					}
				else
					{
					Flat[Base_i_j + HMMSTATE_M] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_IY] = LOG_ZERO;
					Flat[Base_i_j + HMMSTATE_JY] = LOG_ZERO;
					}
				}

			Base_i_j -= BaseInc_j;
			Base_i1_j -= BaseInc_j;
			Base_i_j1 -= BaseInc_j;
			Base_i1_j1 -= BaseInc_j;
			}
		}
	}
