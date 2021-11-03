#include "muscle.h"

float CalcTotalProbFlat(const float *FlatFwd, const float *FlatBwd,
  uint LX, uint LY)
	{
	float Sum = LOG_ZERO;
	uint Base = HMMSTATE_COUNT*(LX*(LY+1) + LY);

	for (uint s = 0; s < HMMSTATE_COUNT; ++s)
		{
		float FwdScore = FlatFwd[Base + s];
		float BwdScore = FlatBwd[Base + s];
		LOG_PLUS_EQUALS(Sum, FwdScore + BwdScore);
		}
	return Sum;
	}
