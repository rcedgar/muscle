#include "muscle.h"
#include "mega.h"

void cmd_test()
	{
	Mega M;
	M.FromFile(g_Arg1);
	const uint FeatureCount = M.GetFeatureCount();
	for (uint i = 0; i < FeatureCount; ++i)
		M.LogFeatureParams(i);
	}
