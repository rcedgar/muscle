#include "muscle.h"
#include "masm.h"

void cmd_masm_train()
	{
	const string &AlnFN = g_Arg1;
	const string &MegaFN = opt(input);

	Mega::FromFile(MegaFN);

	MultiSequence Aln;
	Aln.FromFASTA(AlnFN);

	float GapOpen = 4;
	float GapExt = 0.5;

	MASM M;
	M.FromMSA(Aln, GapOpen, GapExt);
	M.ToFile(opt(output));
	}
