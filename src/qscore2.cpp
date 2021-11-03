#include "muscle.h"
#include "qscorer.h"

void cmd_qscore2()
	{
	const string TestFileName = opt(qscore2);
	const string RefFileName = opt(ref);

	double MaxGapFract = optd(max_gap_fract, 1.0);

	string Name;
	GetBaseName(TestFileName.c_str(), Name);

	MSA Test;
	MSA Ref;
	Test.FromFASTAFile(TestFileName);
	Ref.FromFASTAFile_PreserveCase(RefFileName);

	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;
	QS.Run(Test, Ref);
	ProgressLog("%s: Q=%.4f, TC=%.4f\n", Name.c_str(), QS.m_Q, QS.m_TC);
	}

void cmd_qscoredir()
	{
	const string NamesFileName = opt(qscoredir);
	string TestDir = opt(testdir);
	string RefDir = opt(refdir);
	const string OutputFileName = opt(output);
	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	Dirize(TestDir);
	Dirize(RefDir);

	vector<string> Names;
	ReadStringsFromFile(NamesFileName, Names);

	FILE *fOut = CreateStdioFile(OutputFileName);

	float SumQ = 0;
	float SumTC = 0;

	float AvgQ = 0;
	float AvgTC = 0;

	const uint NameCount = SIZE(Names);
	for (uint i = 0; i < NameCount; ++i)
		{
		ProgressStep(i, NameCount, "%s  Q %.2f TC %.2f",
		  TestDir.c_str(), AvgQ, AvgTC);

		const string &Name = Names[i];
		const string &TestFileName = TestDir + Name;
		const string &RefFileName = RefDir + Name;

		MSA Test;
		MSA Ref;
		Test.FromFASTAFile(TestFileName);

		extern bool g_FASTA_Upper;
		bool SaveUpper = g_FASTA_Upper;
		g_FASTA_Upper = false;
		Ref.FromFASTAFile(RefFileName);
		g_FASTA_Upper = SaveUpper;

		QScorer QS;
		QS.m_MaxGapFract = MaxGapFract;
		QS.Run(Test, Ref);
		Pf(fOut, "set=%s	q=%.4f	tc=%.4f\n", Name.c_str(), QS.m_Q, QS.m_TC); 

		SumQ += QS.m_Q;
		SumTC += QS.m_TC;

		AvgQ = SumQ/(i+1);
		AvgTC = SumTC/(i+1);
		}

	Pf(fOut, "testdir=%s n=%u avgq=%.4f	avgtc=%.4f\n", 
	  TestDir.c_str(), NameCount, AvgQ, AvgTC);
	CloseStdioFile(fOut);
	}
