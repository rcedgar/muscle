#include "muscle.h"
#include "textfile.h"
#include "tree.h"
#include "pprog.h"
#include <map>

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strings);

void PProg::RunGuideTree(const Tree &GuideTree)
	{
	asserta(m_InputMSACount > 0);
	m_JoinCount = m_InputMSACount - 1;
	m_NodeCount = m_InputMSACount + m_JoinCount;

	vector<uint> Indexes1;
	vector<uint> Indexes2;
	GetGuideTreeJoinOrder(GuideTree, m_MSALabelToIndex,
	  Indexes1, Indexes2);

	Run2(Indexes1, Indexes2);
	}

void cmd_pprogt()
	{
	vector<string> MSAFileNames;
	ReadStringsFromFile(opt(pprogt), MSAFileNames);

	const uint MSACount = SIZE(MSAFileNames);
	asserta(MSACount > 1);

	const string &OutputFileName = opt(output);

	PProg PP;
	PP.m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	if (optset_paircount)
		PP.m_TargetPairCount = int(opt(paircount));

	bool IsNucleo;
	PP.LoadMSAs(MSAFileNames, IsNucleo);
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);
	InitProbcons();

	Tree T;
	T.FromFile(opt(guidetreein));

	PP.RunGuideTree(T);

	const MultiSequence &FinalMSA = PP.GetFinalMSA();
	FinalMSA.WriteMFA(OutputFileName);
	PP.WriteGuideTree(opt(guidetreeout));
	}
