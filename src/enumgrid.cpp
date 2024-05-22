#include "myutils.h"

bool GetNextEnumGrid(const vector<uint> &Sizes, vector<uint> &Indexes)
	{
	const uint N = SIZE(Sizes);
	asserta(SIZE(Indexes) == N);
	for (uint i = 0; i < N; ++i)
		{
		uint Index = Indexes[i];
		uint Size = Sizes[i];
		if (Index + 1 < Size)
			{
			Indexes[i] += 1;
			return true;
			}
		Indexes[i] = 0;
		}
	return false;
	}
