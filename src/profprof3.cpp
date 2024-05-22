#include "muscle.h"
#include "profile3.h"
#include "cachemem3.h"
#include "m3alnparams.h"

void cmd_profprof3()
	{
	const string &InputFileName = opt(profprof3);
	const string &InputFileName2 = opt(input2);
	M3AlnParams AP;
	AP.SetFromCmdLine(false);
	const Mx2020 &SubstMx_Letter =
	  AP.m_SubstMx_Letter;
	float GapOpen = AP.m_GapOpen;

	MultiSequence MSA1;
	MultiSequence MSA2;
	MSA1.FromFASTA(InputFileName);
	MSA2.FromFASTA(InputFileName2);
	bool IsNucleo = MSA1.GuessIsNucleo();
	bool IsNucleo2 = MSA2.GuessIsNucleo();
	asserta(IsNucleo == IsNucleo2);
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);

	const uint SeqCount1 = MSA1.GetSeqCount();
	const uint SeqCount2 = MSA2.GetSeqCount();
	float w1 = 1.0f/SeqCount1;
	float w2 = 1.0f/SeqCount2;
	vector<float> SeqWeights1(SeqCount1, w1);
	vector<float> SeqWeights2(SeqCount2, w2);

	const uint SeqCount12 = SeqCount1 + SeqCount2;
	float w12 = 1.0f/SeqCount12;
	vector<float> SeqWeights12(SeqCount12, w12);

	Profile3 Prof1;
	Profile3 Prof2;
	Prof1.FromMSA(MSA1, SubstMx_Letter, GapOpen, SeqWeights1);
	Prof2.FromMSA(MSA2, SubstMx_Letter, GapOpen, SeqWeights2);
	Log("_____________ Prof1 ________________________\n");
	Prof1.LogMe(&MSA1);
	Log("_____________ Prof2 ________________________\n");
	Prof2.LogMe(&MSA2);
	Prof1.Validate();
	Prof2.Validate();
	Prof1.ToTsv(opt(output1));
	Prof2.ToTsv(opt(output2));

	CacheMem3 CM;
	string Path;
	float Score = NWSmall3(CM, Prof1, Prof2, Path);
	Log("Score=%.4g\n", Score);
	Log("Path=%s\n", Path.c_str());

	MultiSequence MSA12;
	AlignTwoMSAsGivenPath(MSA1, MSA2, Path, MSA12);
	MSA12.WriteMFA(opt(output));

	Profile3 Prof12Msa;
	Prof12Msa.FromMSA(MSA12, SubstMx_Letter, GapOpen, SeqWeights12);

	Log("\n______________________ Prof12Msa _________________________\n");
	Prof12Msa.LogMe(&MSA12);
	Prof12Msa.ToTsv(opt(output4));
	Prof12Msa.Validate();

	Profile3 Prof12Path;
	AlignTwoProfsGivenPath(Prof1, 0.5f, Prof2, 0.5f,
	  SubstMx_Letter, GapOpen, Path, Prof12Path);
	Prof12Path.ToTsv(opt(output3));
	Log("\n______________________ Prof12Path _________________________\n");
	Prof12Path.LogMe(&MSA12);
	Prof12Path.Validate();

	uint DiffCount = Prof12Msa.LogDiffs(Prof12Path);
	ProgressLog("%u diffs\n", DiffCount);

	Prof12Path.ToTsv("prof5.tsv");
	}
