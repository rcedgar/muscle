#include "muscle.h"
#include "qscorer2.h"

void cmd_qscore2()
	{
	const string TestFileName = opt(qscore2);
	const string RefFileName = opt(ref);

	double MaxGapFract = optd(max_gap_fract, 1.0);

	string Name;
	GetBaseName(TestFileName.c_str(), Name);

	MSA Test;
	Test.FromFASTAFile(TestFileName);

	MultiSequence Ref;
	Ref.m_DupeLabelsOk = true;
	Ref.FromFASTA(RefFileName);

	QScorer2 QS;
	double Q = QS.Run(Test, Ref);

	ProgressLog("Q=%.4f %s\n", Q, Name.c_str());
	}

void cmd_qscore2dir()
	{
	const string NamesFileName = g_Arg1;
	string TestDir = opt(testdir);
	string RefDir = opt(refdir);
	const string OutputFileName = opt(output);

	string Algo;
	Algo = (string) BaseName(TestDir.c_str());

	Dirize(TestDir);
	Dirize(RefDir);

	vector<string> Names;
	ReadStringsFromFile(NamesFileName, Names);

	FILE *fOut = CreateStdioFile(OutputFileName);

	double SumQ = 0;
	double AvgQ = 0;

	const uint NameCount = SIZE(Names);
	for (uint i = 0; i < NameCount; ++i)
		{
		ProgressStep(i, NameCount, "%s  Q %.3f",
		  TestDir.c_str(), AvgQ);

		const string &Name = Names[i];
		const string &TestFileName = TestDir + Name;
		const string &RefFileName = RefDir + Name;

		MSA Test;
		MultiSequence Ref;
		Test.FromFASTAFile(TestFileName);

		extern bool g_FASTA_Upper;
		bool SaveUpper = g_FASTA_Upper;
		g_FASTA_Upper = false;
		Ref.m_DupeLabelsOk = true;
		Ref.FromFASTA(RefFileName);
		g_FASTA_Upper = SaveUpper;

		QScorer2 QS;
		double Q = QS.Run(Test, Ref);
		Pf(fOut, "set=%s	q=%.4f\n", Name.c_str(), Q); 

		SumQ += Q;
		AvgQ = SumQ/(i+1);
		}

	Pf(fOut, "testdir=%s n=%u avgq=%.4f\n", 
	  TestDir.c_str(), NameCount, AvgQ);
	Pf(fOut, "Algo=%s Q=%.3f\n", Algo.c_str(), AvgQ);
	CloseStdioFile(fOut);

	printf("Algo=%s Q=%.3f\n", Algo.c_str(), AvgQ);
	}
