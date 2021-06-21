#include "muscle.h"
#include "qscorer.h"

void cmd_qscore2()
	{
	const string TestFileName = opt(qscore2);
	const string RefFileName = opt(ref);
	const string ReportFileName = opt(report);

	string Name;
	GetBaseName(TestFileName.c_str(), Name);

	MSA Test;
	MSA Ref;
	Test.FromFASTAFile(TestFileName);
	Ref.FromFASTAFile_PreserveCase(RefFileName);

	FILE *fRep = CreateStdioFile(ReportFileName);

	QScorer QS;
	QS.Run(TestFileName, RefFileName, Test, Ref);
	QS.Report(fRep);
	ProgressLog("%s: Q=%.4f, TC=%.4f, CF=%.3g\n", Name.c_str(), QS.m_Q, QS.m_TC, QS.m_CF);

	CloseStdioFile(fRep);
	}

void cmd_qscoredir()
	{
	const string NamesFileName = opt(qscoredir);
	string TestDir = opt(testdir);
	string RefDir = opt(refdir);
	const string OutputFileName = opt(output);

	Dirize(TestDir);
	Dirize(RefDir);

	vector<string> Names;
	ReadStringsFromFile(NamesFileName, Names);

	FILE *fOut = CreateStdioFile(OutputFileName);

	float SumQ = 0;
	float SumCF = 0;
	float SumTC = 0;
	float SumTC90 = 0;
	float SumAcc = 0;

	float AvgQ = 0;
	float AvgCF = 0;
	float AvgTC = 0;
	float AvgTC90 = 0;
	float AvgAcc = 0;

	const uint NameCount = SIZE(Names);
	for (uint i = 0; i < NameCount; ++i)
		{
		ProgressStep(i, NameCount, "%s  Q %.2f TC %.2f CF %.2f TC90 %.2f Acc %.2f",
		  TestDir.c_str(), AvgQ, AvgTC, AvgCF, AvgTC90, AvgAcc);

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
		QS.Run(TestFileName, RefFileName, Test, Ref);
		Pf(fOut, "set=%s	q=%.4f	tc=%.4f	cf=%.4f	tc90=%.4f	acc=%.4f\n",
		  Name.c_str(), QS.m_Q, QS.m_TC, QS.m_CF, QS.m_TC90, QS.m_Acc); 

		SumQ += QS.m_Q;
		SumTC += QS.m_TC;
		SumCF += QS.m_CF;
		SumTC90 += QS.m_TC90;
		SumAcc += QS.m_Acc;

		AvgQ = SumQ/(i+1);
		AvgTC = SumTC/(i+1);
		AvgTC90 = SumTC90/(i+1);
		AvgCF = SumCF/(i+1);
		AvgAcc = SumAcc/(i+1);
		}

	Pf(fOut, "testdir=%s n=%u avgq=%.4f	avgtc=%.4f	avgcf=%.4f	avgtc90=%.4f	avgacc=%.4f\n", 
	  TestDir.c_str(), NameCount, AvgQ, AvgTC, AvgCF, AvgTC90, AvgAcc);
	CloseStdioFile(fOut);
	}
