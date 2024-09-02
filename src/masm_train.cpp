#include "muscle.h"
#include "masm.h"

void cmd_masm_stats()
	{
	MASM M;
	M.FromFile(g_Arg1);
	ProgressLog("%10u  Sequences\n", M.m_SeqCount);
	ProgressLog("%10u  Columns\n", M.m_ColCount);
	ProgressLog("%10u  Features ", M.m_FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < M.m_FeatureCount; ++ FeatureIdx)
		ProgressLog(" %s/%u",
		  M.m_FeatureNames[FeatureIdx].c_str(), 
		  M.m_AlphaSizes[FeatureIdx]); 
	ProgressLog("\n");
	}

void cmd_masm_train()
	{
	const string &AlnFN = g_Arg1;
	const string &MegaFN = opt(input);

	Mega::FromFile(MegaFN);

	MultiSequence Aln;
	Aln.FromFASTA(AlnFN);

	string Label;
	if (optset_label)
		Label = opt(label);
	else
		Label = string(BaseName(AlnFN.c_str()));

	MASM M;
	M.FromMSA(Aln, Label, Mega::m_GapOpen, Mega::m_GapExt);
	M.ToFile(opt(output));
	}
