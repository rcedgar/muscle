#include "muscle.h"

void LogDistMx(const string &Msg, const vector<vector<float> > &Mx)
	{
	Log("\n");
	Log("LogDistMx(%s)\n", Msg.c_str());
	const uint RowCount = SIZE(Mx);
	asserta(RowCount > 0);
	const uint ColCount = SIZE(Mx[0]);
	for (uint Row = 0; Row < RowCount; ++Row)
		{
		Log("[%5u]  ", Row);
		const uint ColCount = SIZE(Mx[Row]);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			float x = Mx[Row][Col];
			if (x == FLT_MAX)
				Log("  %7.7s", "*");
			else if (x == LOG_ZERO)
				Log("  %7.7s", ".");
			else
				Log("  %7.3g", x);
			}
		Log("\n");
		}
	}
