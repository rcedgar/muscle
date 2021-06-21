#include "myutils.h"
#include "muscle.h"
#include "probcons.h"
#include "pprog.h"

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strings);

void PProg::AlignAndJoin(uint Index1, uint Index2)
	{
	m_JoinMSAIndexes1.push_back(Index1);
	m_JoinMSAIndexes2.push_back(Index2);

	MultiSequence &MSA1 = (MultiSequence &) GetMSA(Index1);
	MultiSequence &MSA2 = (MultiSequence &) GetMSA(Index2);
	AssertSameLabels(MSA1);
	AssertSameLabels(MSA2);

	const string &MSALabel1 = m_MSALabels[Index1];
	const string &MSALabel2 = m_MSALabels[Index2];
	string ProgressStr;
	Ps(ProgressStr, "Join %u / %u", m_JoinIndex+1, m_JoinCount);

	vector<char> Path;
	float Score = AlignMSAs(ProgressStr, MSA1, MSA2, m_TargetPairCount, Path);

	string MSALabel12;
	Ps(MSALabel12, "Join_%u", m_JoinIndex+1);

	MultiSequence *MSA12 = new MultiSequence;
	AlignMSAsByPath(MSA1, MSA2, Path, *MSA12);
	AssertSeqsEq(MSA1, *MSA12);
	AssertSeqsEq(MSA2, *MSA12);
	AssertSameSeqsJoin(MSA1, MSA2, *MSA12);
	AssertSameLabels(*MSA12);

	uint NewMSAIndex = m_InputMSACount + m_JoinIndex;
	SetMSA(NewMSAIndex, *MSA12);
	if (optset_savedir)
		{
		string Prefix = opt(savedir);
		Dirize(Prefix);
		string JoinFileName;
		Ps(JoinFileName, "%sjoin%u", Prefix.c_str(), m_JoinIndex);
		ProgressLog("Writing join MSA: %s\n", JoinFileName.c_str());
		MSA12->WriteMFA(JoinFileName);
		}
	}

void PProg::Run2(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2)
	{
	asserta(m_InputMSACount > 0);
	m_JoinCount = m_InputMSACount - 1;
	m_NodeCount = m_InputMSACount + m_JoinCount;

	asserta(SIZE(Indexes1) == m_JoinCount);
	asserta(SIZE(Indexes2) == m_JoinCount);

	ValidateJoinOrder(Indexes1, Indexes2);

	for (m_JoinIndex = 0; m_JoinIndex < m_JoinCount; ++m_JoinIndex)
		{
		uint Index1 = Indexes1[m_JoinIndex];
		uint Index2 = Indexes2[m_JoinIndex];
		AlignAndJoin(Index1, Index2);
		}
	}

void cmd_pprog2()
	{
	PProg PP;
	vector<string> MSAFileNames;
	ReadStringsFromFile(opt(pprog2), MSAFileNames);

	const uint MSACount = SIZE(MSAFileNames);
	asserta(MSACount > 1);
	const string &OutputFileName = opt(output);

	PP.m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	if (optset_paircount)
		PP.m_TargetPairCount = int(opt(paircount));

	InitProbcons();

	vector<uint> Indexes1;
	vector<uint> Indexes2;
	FILE *f = OpenStdioFile(opt(joins));
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		uint Index1 = StrToUint(Fields[0]);
		uint Index2 = StrToUint(Fields[1]);
		Indexes1.push_back(Index1);
		Indexes2.push_back(Index2);
		}
	CloseStdioFile(f);
	f = 0;

	PP.LoadMSAs(MSAFileNames);
	PP.Run2(Indexes1, Indexes2);

	const MultiSequence &FinalMSA = PP.GetFinalMSA();
	FinalMSA.WriteMFA(OutputFileName);
	}
