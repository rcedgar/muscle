#pragma once

static void FlatCoords(uint Ix, uint LY, uint &s, uint &i, uint &j)
	{
	s = Ix%HMMSTATE_COUNT;
	uint r = (Ix - s)/HMMSTATE_COUNT;	 // i*(LY+1) + j
	j = r%(LY+1);
	i = r/(LY+1);
	}

static inline uint FlatIx(uint s, uint i, uint j, uint LY)
	{
	uint Ix = HMMSTATE_COUNT*(i*(LY+1) + j) + s;
	return Ix;
	}

static inline uint FlatIx(HMMSTATE s, uint i, uint j, uint LY)
	{
	uint Ix = FlatIx(uint(s), i, j, LY);
	return Ix;
	}
