#include "myutils.h"
#include "muscle.h"
#include "textfile.h"
#include "tree.h"
#include "pprog.h"
#include <map>

void cmd_pprog_tree()
	{
	//const string &InputFileName = opt(pprog_tree);
	asserta(optset_guidetreein);
	const string &OutputFileName = opt(output);

	MultiSequence InputSeqs;
	LoadInput(InputSeqs);
	//InputSeqs.FromFASTA(InputFileName);
	bool IsNucleo = InputSeqs.GuessIsNucleo();

	const uint SeqCount = InputSeqs.GetSeqCount();
	vector<const MultiSequence *> MSAs;
	vector<string> MSALabels;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = InputSeqs.GetSequence(i);
		const string &Label = InputSeqs.GetLabel(i);
		MultiSequence *MS = new MultiSequence;
		MS->AddSequence(Seq, false);
		MSAs.push_back(MS);
		MSALabels.push_back(Label);
		}

	PProg PP;
	PP.m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	if (optset_paircount)
		PP.m_TargetPairCount = int(opt(paircount));

	PP.SetMSAs(MSAs, MSALabels);
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);
	InitProbcons();

	Tree T;
	T.FromFile(opt(guidetreein));

	PP.RunGuideTree(T);

	const MultiSequence &FinalMSA = PP.GetFinalMSA();
	FinalMSA.WriteMFA(OutputFileName);
	}
