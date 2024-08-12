#include "muscle.h"
#include "kmerdist33.h"
#include "kmerdist66.h"
#include "clustalweights.h"
#include "pprog3.h"
#include "muscle3.h"

void cmd_muscle3()
	{
	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(g_Arg1);
	bool IsNucleo = InputSeqs.GuessIsNucleo();

	M3AlnParams AP;
	AP.SetFromCmdLine(IsNucleo);
	Muscle3 M3;
	M3.m_AP = &AP;
	if (opt(muscle3_randomorder))
		M3.RunRO(AP, InputSeqs);
	else
		{
		M3.Run(AP, InputSeqs);
		if (optset_guidetreeout)
			M3.m_GuideTree.ToFile(opt(guidetreeout));
		}
	M3.WriteMSA(opt(output));
	}
