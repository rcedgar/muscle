#include "myutils.h"
#include "diagbox.h"

#define TEST	0

/***
DiagBox represents a diagonal "rectangle" in the D.P. matrix.

	i = 0..LA-1
	j = 0..LB-1
	d = LA - i + j = 1 .. LA+LB-1
	j = d - LA + i
	i = LA - d + j
***/

void GetDiagRange(uint LA, uint LB, uint d,
  uint &mini, uint &minj, uint &maxi, uint &maxj)
	{
	if (d >= LA)
		{
		mini = 0;
		maxi = min(LA+LB-1-d, LA-1);
		minj = d - LA;
		maxj = min(LB-1, d-1);
		}
	else
		{
		mini = LA-d;
		maxi = min(LA+LB-1-d, LA-1);
		minj = 0;
		maxj = min(LB-1, d-1);
		}
	}

void GetDiagBox(uint LA, uint LB, uint DiagLo, uint DiagHi, DiagBox &Box)
	{
	asserta(DiagLo <= DiagHi);
	asserta(DiagLo >= 1);
	asserta(DiagHi <= LA + LB - 1);

	Box.LA = LA;
	Box.LB = LB;

	Box.dlo = DiagLo;
	Box.dhi = DiagHi;

	GetDiagRange(LA, LB, DiagLo, Box.dlo_mini, Box.dlo_minj, Box.dlo_maxi, Box.dlo_maxj);
	GetDiagRange(LA, LB, DiagHi, Box.dhi_mini, Box.dhi_minj, Box.dhi_maxi, Box.dhi_maxj);
	}

void GetDiagLoHi(uint LA, uint LB, const char *Path,
  uint &dlo, uint &dhi)
	{
	dlo = UINT_MAX;
	dhi = UINT_MAX;

	uint i = 0;
	uint j = 0;
	for (uint k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			break;
		if (c == 'M')
			{
			uint d = LA - i + j;
			if (dlo == UINT_MAX)
				{
				dlo = d;
				dhi = d;
				}
			else
				{
				if (d < dlo)
					dlo = d;
				if (d > dhi)
					dhi = d;
				}
			}
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	}

#if TEST
static void Test2(uint LA, uint LB, uint DiagLo, uint DiagHi)
	{
	DiagBox Box;
	GetDiagBox(LA, LB, DiagLo, DiagHi, Box);
	Box.LogMe();
	Box.Validate();
	}

static void Test1(uint LA, uint LB, uint d,
  uint i, uint j, uint I, uint J)
	{
	uint mini, maxi, minj, maxj;
	GetDiagRange(LA, LB, d, mini, minj, maxi, maxj);
	Log("LA=%u LB=%u d=%u (%u,%u) (%u,%u) expected (%u,%u) (%u,%u)\n",
	  LA, LB, d, mini, minj, maxi, maxj, i, j, I, J);

	asserta(mini == i);
	asserta(maxi == I);
	asserta(minj == j);
	asserta(maxj == J);
	}

void TestDiagBox()
	{
	Test2(16, 19, 17, 37);

	Test1(5, 3, 1, 4, 0, 4, 0);
	Test1(5, 3, 2, 3, 0, 4, 1);
	Test1(5, 3, 3, 2, 0, 4, 2);
	Test1(5, 3, 4, 1, 0, 3, 2);
	Test1(5, 3, 5, 0, 0, 2, 2);
	Test1(5, 3, 6, 0, 1, 1, 2);
	Test1(5, 3, 7, 0, 2, 0, 2);

	Test1(3, 5, 1, 2, 0, 2, 0);
	Test1(3, 5, 2, 1, 0, 2, 1);
	Test1(3, 5, 3, 0, 0, 2, 2);
	Test1(3, 5, 4, 0, 1, 2, 3);
	Test1(3, 5, 5, 0, 2, 2, 4);
	Test1(3, 5, 6, 0, 3, 1, 4);
	Test1(3, 5, 7, 0, 4, 0, 4);

	Test1(5, 5, 1, 4, 0, 4, 0);
	Test1(5, 5, 2, 3, 0, 4, 1);
	Test1(5, 5, 3, 2, 0, 4, 2);
	Test1(5, 5, 4, 1, 0, 4, 3);
	Test1(5, 5, 5, 0, 0, 4, 4);
	Test1(5, 5, 6, 0, 1, 3, 4);
	Test1(5, 5, 7, 0, 2, 2, 4);
	Test1(5, 5, 8, 0, 3, 1, 4);
	Test1(5, 5, 9, 0, 4, 0, 4);

	for (uint LA = 2; LA <= 5; ++LA)
		for (uint LB = 2; LB <= 5; ++LB)
			for (uint dlo = 1; dlo <= LA+LB-1; ++dlo)
				for (uint dhi = dlo; dhi <= LA+LB-1; ++dhi)
					Test2(LA, LB, dlo, dhi);

	Log("\n");
	Log("ALL OK\n");
	}
#endif // TEST
