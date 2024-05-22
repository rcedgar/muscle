#include "myutils.h"

void GetRandomDistMx(uint SeqCount, vector<vector<float> > &DistMx)
	{
	DistMx.clear();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount);
		DistMx[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			float d = (randu32()%100)/100.0f + 0.001f;
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}
