#include "muscle.h"

void TraceBackFlat(const char *TB, uint LX, uint LY, string &Path)
	{
	Path.clear();
	int i = int(LX);
	int j = int(LY);
	for (;;)
		{
		if (i == 0 && j == 0)
			break;
		if (i < 0 || j < 0)
			{
			Warning("TraceBackFlat i=%d j=%d", i, j);
			return;
			}
		char TBChar = TB[i*(LY+1) + j];
		Path.push_back(TBChar);

		switch (TBChar)
			{
		case 'B':
			--i;
			--j;
			break;

		case 'X':
			--i;
			break;

		case 'Y':
			--j;
			break;
			}
		}
	reverse(Path.begin(), Path.end());
	}
