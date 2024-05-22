#include "muscle.h"
#include "qscorer.h"

static void GetCommonLabels(const MultiSequence &Test1, 
  const MultiSequence &Test2, set<string> &CommonLabels)
	{
	CommonLabels.clear();
	const uint SeqCount1 = Test1.GetSeqCount();
	const uint SeqCount2 = Test2.GetSeqCount();

	set<string> Labels1;
	for (uint i = 0; i < SeqCount1; ++i)
		Labels1.insert((string) Test1.GetLabel(i));

	for (uint i = 0; i < SeqCount2; ++i)
		{
		const string &Label2 = (string) Test2.GetLabel(i);
		if (Labels1.find(Label2) != Labels1.end())
			CommonLabels.insert(Label2);
		}
	}

static void MakeRefCommon(const MultiSequence &Ref, const set<string> &CommonLabels,
  MultiSequence  &RefCommon)
	{
	RefCommon.Clear();
	map<string, uint> LabelToIndex;
	for (uint i = 0; i < Ref.GetSeqCount(); ++i)
		{
		string Label = (string) Ref.GetLabel(i);
		LabelToIndex[Label] = i;
		}

	for (set<string>::const_iterator iter = CommonLabels.begin();
	  iter != CommonLabels.end(); ++iter)
		{
		const string &Label = *iter;
		map<string, uint>::const_iterator iter2 = LabelToIndex.find(Label);
		if (iter2 == LabelToIndex.end())
			continue;
		uint SeqIndex = iter2->second;
		const Sequence *Seq = Ref.GetSequence(SeqIndex);
		RefCommon.AddSequence(Seq, false);
		}
	}

void cmd_qscoredir_cmp2()
	{
	const string NamesFileName = g_Arg1;
	string TestDir1 = opt(testdir);
	string TestDir2 = opt(input2);
	string RefDir = opt(refdir);
	const string OutputFileName = opt(output);

	string Algo1;
	string Algo2;
	Algo1 = (string) BaseName(TestDir1.c_str());
	Algo2 = (string) BaseName(TestDir2.c_str());

	Dirize(TestDir1);
	Dirize(TestDir2);
	Dirize(RefDir);

	string Name1;
	string Name2;
	GetBaseName(TestDir1, Name1);
	GetBaseName(TestDir2, Name2);

	vector<string> Names;
	ReadStringsFromFile(NamesFileName, Names);

	FILE *fOut = CreateStdioFile(OutputFileName);

	double SumQ1 = 0;
	double SumQ2 = 0;
	double SumTC1 = 0;
	double SumTC2 = 0;

	double AvgQ1 = 0;
	double AvgQ2 = 0;
	double AvgTC1 = 0;
	double AvgTC2 = 0;

	const uint NameCount = SIZE(Names);
	uint N = 0;
	for (uint i = 0; i < NameCount; ++i)
		{
		const string &Name = Names[i];
		ProgressStep(i, NameCount, "%s  Q %.2f,%.2f TC %.2f,%.2f",
		  Name.c_str(), AvgQ1, AvgQ2, AvgTC1, AvgTC2);

		const string &TestFileName1 = TestDir1 + Name;
		const string &TestFileName2 = TestDir2 + Name;
		const string &RefFileName = RefDir + Name;

		bool Exists1 = StdioFileExists(TestFileName1);
		bool Exists2 = StdioFileExists(TestFileName2);
		if (!Exists1)
			Log("Missing file1 %s\n", TestFileName1.c_str());
		if (!Exists2)
			Log("Missing file2 %s\n", TestFileName2.c_str());
		if (!Exists1 || !Exists2)
			continue;

		MultiSequence Test1;
		MultiSequence Test2;
		MultiSequence Ref;
		Test1.FromFASTA(TestFileName1);
		Test2.FromFASTA(TestFileName2);

		const uint n1 = Test1.GetSeqCount();
		const uint n2 = Test2.GetSeqCount();

		extern bool g_FASTA_Upper;
		bool SaveUpper = g_FASTA_Upper;
		g_FASTA_Upper = false;
		Ref.FromFASTA(RefFileName);
		g_FASTA_Upper = SaveUpper;
		const uint n = Ref.GetSeqCount();

		set<string> CommonLabels;
		GetCommonLabels(Test1, Test2, CommonLabels);
		if (CommonLabels.size() < 2)
			continue;

		MultiSequence RefCommon;
		MakeRefCommon(Ref, CommonLabels, RefCommon);
		if (RefCommon.GetSeqCount() < 2)
			continue;

		QScorer QS;
		QS.m_MissingTestSeqOk = true;
		QS.Run(Name, Test1, RefCommon);
		double Q1 = QS.m_Q;
		double TC1 = QS.m_TC;

		QS.Run(Name, Test2, RefCommon);
		double Q2 = QS.m_Q;
		double TC2 = QS.m_TC;

		Pf(fOut, "set=%s\tq1=%.4f\tq2=%.4f\ttc1=%.4f\ttc2=%.4f\tn=%u\tn1=%u\tn2=%u\n",
		  Name.c_str(), Q1, Q2, TC1, TC2, n, n1, n2);

		asserta(!isnan(QS.m_Q));
		if (QS.m_Q >= 0)
			{
			++N;
			SumQ1 += Q1;
			SumQ2 += Q2;
			SumTC1 += TC1;
			SumTC2 += TC2;

			AvgQ1 = SumQ1/N;
			AvgQ2 = SumQ2/N;
			AvgTC1 = SumTC1/N;
			AvgTC2 = SumTC2/N;
			}
		}

	Pf(fOut, "testdir1=%s\ttestdir2=%s\tn=%u\tN=%d\tavgq1=%.4f\tavgq2=%.4f\tavgtc1=%.4f\tavgtc2=%.4f\n", 
	  TestDir1.c_str(), TestDir2.c_str(),
	  N, NameCount, AvgQ1, AvgQ2, AvgTC1, AvgTC2);
	CloseStdioFile(fOut);
	}
