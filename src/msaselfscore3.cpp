#include "muscle.h"
#include "locallock.h"
#include "m3alnparams.h"

void cmd_msaselfscore3()
	{
	const string &InputFileName = opt(msaselfscore3);

	MultiSequence MSA;
	MSA.FromFASTA(InputFileName);
	bool IsNucleo = MSA.GuessIsNucleo();
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);

	M3AlnParams AP;
	AP.SetFromCmdLine(IsNucleo);
	Mx2020 &SubstMx_Letter = AP.m_SubstMx_Letter;
	float GapOpen = AP.m_GapOpen;

	const uint SeqCount = MSA.GetSeqCount();
	float w = 1.0f/SeqCount;
	vector<float> SeqWeights(SeqCount, w);

	Profile3 Prof;
	Prof.FromMSA(MSA, SubstMx_Letter, GapOpen, SeqWeights);
	Prof.Validate();
	float SelfScore = Prof.GetSelfScore();

	ProgressLog("MSASelfScore3=%.5g, MSA=%s\n", 
	  SelfScore, BaseName(InputFileName.c_str()));
	}
