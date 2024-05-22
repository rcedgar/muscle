#include "muscle.h"
#include "qscorer.h"

void cmd_qscore()
	{
	const string TestFileName = opt(qscore);
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
	QS.Run(Name, Test, Ref);
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

	string Algo;
	Algo = (string) BaseName(TestDir.c_str());

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
	uint N = 0;
	uint MissingTestFileCount = 0;
	for (uint i = 0; i < NameCount; ++i)
		{
		ProgressStep(i, NameCount, "%s  Q %.2f TC %.2f",
		  TestDir.c_str(), AvgQ, AvgTC);

		const string &Name = Names[i];
		const string &TestFileName = TestDir + Name;
		const string &RefFileName = RefDir + Name;

		MSA Test;
		MSA Ref;
		if (opt(missingtestfileok) && !StdioFileExists(TestFileName))
			{
			++MissingTestFileCount;
			Log("Missing test file %s\n", TestFileName.c_str());
			continue;
			}
		Test.FromFASTAFile(TestFileName);

		extern bool g_FASTA_Upper;
		bool SaveUpper = g_FASTA_Upper;
		g_FASTA_Upper = false;
		Ref.FromFASTAFile(RefFileName);
		g_FASTA_Upper = SaveUpper;

		QScorer QS;
		QS.m_MaxGapFract = MaxGapFract;
		QS.Run(Name, Test, Ref);
		Pf(fOut, "set=%s	q=%.4f	tc=%.4f\n", Name.c_str(), QS.m_Q, QS.m_TC); 

		asserta(!isnan(QS.m_Q));
		if (QS.m_Q >= 0)
			{
			++N;
			SumQ += QS.m_Q;
			SumTC += QS.m_TC;

			AvgQ = SumQ/N;
			AvgTC = SumTC/N;
			}
		}

	Pf(fOut, "testdir=%s n=%u avgq=%.4f	avgtc=%.4f\n", 
	  TestDir.c_str(), NameCount, AvgQ, AvgTC);
	Pf(fOut, "Algo=%s Q=%.3f\n", Algo.c_str(), AvgQ);
	CloseStdioFile(fOut);

	printf("Algo=%s Q=%.3f N=%u/%u\n", Algo.c_str(), AvgQ, N, NameCount);

	if (MissingTestFileCount > 0)
		Warning("%u missing test MSAs (see log)", MissingTestFileCount);
	}
