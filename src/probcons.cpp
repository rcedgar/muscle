#include "muscle.h"

void ProgressLogInputSummary(const string &FileName, const MultiSequence &Seqs)
	{
	const uint SeqCount = (uint) Seqs.GetNumSequences();
	uint MinL = 0;
	uint MaxL = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence &Seq = *Seqs.GetSequence(i);
		uint L = (uint) Seq.GetLength();
		if (i == 0)
			{
			MinL = L;
			MaxL = L;
			}
		else
			{
			MinL = min(MinL, L);
			MaxL = max(MaxL, L);
			}
		}
	ProgressLog("\n");
	ProgressLog("Input  %s\n", FileName.c_str());
	ProgressLog("Seqs   %u\n", SeqCount);
	ProgressLog("MinL   %u\n", MinL);
	ProgressLog("MaxL   %u\n", MaxL);
	ProgressLog("\n");
	}

void ProgressLogMSASummary(const string &Str, const MultiSequence &MSA)
	{
	const uint SeqCount = (uint) MSA.GetNumSequences();
	const uint ColCount = (uint) MSA.GetColCount();
	ProgressLog("\n");
	ProgressLog("%s\n", Str.c_str());
	ProgressLog("Seqs   %u\n", SeqCount);
	ProgressLog("Cols   %u\n", ColCount);
	ProgressLog("\n");
	}

//void RunMTProbcons(MultiSequence &InputSeqs)
//	{
//	bool IsNucleo = InputSeqs.GuessIsNucleo();
//	if (IsNucleo)
//		SetAlpha(ALPHA_Nucleo);
//	else
//		SetAlpha(ALPHA_Amino);
//
//	const string &OutputFileName = opt(output);
//
//	InitProbcons();
//
//	MultiSequence* alignment = RunMPC(&InputSeqs);
//
//	alignment->WriteMFA(opt(output));
//
//	string s;
//	Psa(s, "MSA %s", OutputFileName.c_str());
//	ProgressLogMSASummary(s, *alignment);
//	}
