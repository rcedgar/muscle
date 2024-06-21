#include "muscle.h"
#include "mega.h"

void cmd_test()
	{
	Mega M;
	M.FromFile(g_Arg1);
	M.LogFeatureParams(0);
	}
