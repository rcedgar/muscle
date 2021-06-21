#include "muscle.h"
#include "qscorer.h"
#include "qscorer3.h"

void QScorer3::Run(const string &TestFileName1, const string &TestFileName2,
  const string &RefFileName)
	{
	m_Test1.FromFASTAFile(TestFileName1);
	m_Test2.FromFASTAFile(TestFileName2);

	m_Ref.FromFASTAFile_PreserveCase(RefFileName);

	m_QS1.Run(TestFileName1, RefFileName, m_Test1, m_Ref);
	m_QS2.Run(TestFileName2, RefFileName, m_Test2, m_Ref);

	m_Indexes2.clear();

	asserta(m_QS1.m_RefCols == m_QS2.m_RefCols);
	asserta(m_QS1.m_RefAlignedColCount == m_QS2.m_RefAlignedColCount);

	m_RefCols = &m_QS1.m_RefCols;
	m_RefAlignedColCount = m_QS1.m_RefAlignedColCount;

	map<string, uint> LabelToIndex2;
	for (uint i = 0; i < SIZE(m_QS2.m_Labels); ++i)
		{
		const string &Label2 = m_QS2.m_Labels[i];
		LabelToIndex2[Label2] = i;
		}

	for (uint Index = 0; Index < SIZE(m_QS1.m_Labels); ++Index)
		{
		const string &Label = m_QS1.m_Labels[Index];
		uint RefSeqIndex1 = m_QS1.m_RefSeqIndexes[Index];
		const string &RefLabel = (const string &) m_QS1.m_Ref->GetSeqName(RefSeqIndex1);
		asserta(Label == RefLabel);

		map<string, uint>::const_iterator p =
		  LabelToIndex2.find(Label);
		asserta(p != LabelToIndex2.end());
		uint Index2 = p->second;
		m_Indexes2.push_back(Index2);
		}

	TransQ();
	}

void QScorer3::ToTSV(const string &FileName) const
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	ToTSV(f);
	CloseStdioFile(f);
	}

void QScorer3::ToTSV(FILE *f) const
	{
	if (f == 0)
		return;
	const uint PairCount = SIZE(m_Pairs);
	asserta(SIZE(m_PairIndexToQ1) == PairCount);
	asserta(SIZE(m_PairIndexToQ2) == PairCount);
	asserta(SIZE(m_PairIndexToPWC) == PairCount);
	float SumQ1 = 0;
	float SumQ2 = 0;
	float SumPWC = 0;
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		const pair<uint, uint> &Pair = m_Pairs[PairIndex];
		uint Indexi = Pair.first;
		uint Indexj = Pair.second;
		const string &Labeli = m_QS1.m_Labels[Indexi];
		const string &Labelj = m_QS1.m_Labels[Indexj];
		float Q1 = m_PairIndexToQ1[PairIndex];
		float Q2 = m_PairIndexToQ2[PairIndex];
		float PWC = m_PairIndexToPWC[PairIndex];
		SumQ1 += Q1;
		SumQ2 += Q2;
		SumPWC += PWC;
		if (!opt(allpairs))
			continue;

		fprintf(f, "@PAIR");
		fprintf(f, "	Q1=%.4f", Q1);
		fprintf(f, "	Q2=%.4f", Q2);
		fprintf(f, "	PWC=%.4f", PWC);
		fprintf(f, "	%s",	Labeli.c_str());
		fprintf(f, "	%s",	Labelj.c_str());
		fprintf(f, "	%s",	BaseName(m_QS1.m_RefName.c_str()));
		fprintf(f, "\n");
		}

	float AvgQ1 = SumQ1/PairCount;
	float AvgQ2 = SumQ2/PairCount;
	float AvgPWC = SumPWC/PairCount;

	fprintf(f, "@SET");
	fprintf(f, "	Q1=%.4f", AvgQ1);
	fprintf(f, "	Q2=%.4f", AvgQ2);
	fprintf(f, "	PWC=%.4f", AvgPWC);
	fprintf(f, "	%s",	BaseName(m_QS1.m_RefName.c_str()));
	fprintf(f, "\n");
	}

void cmd_qscore3()
	{
	const string TestFileName1 = opt(qscore3);
	const string TestFileName2 = opt(input2);
	const string RefFileName = opt(ref);
	const string OutputFileName = opt(output);

	QScorer3 Q3;
	Q3.Run(TestFileName1, TestFileName2, RefFileName);
	Q3.ToTSV(OutputFileName);
	}

void cmd_qscore3dir()
	{
	const string SetsFileName = opt(qscore3dir);
	string &TestDir1 = opt(testdir1);
	string &TestDir2 = opt(testdir2);
	string RefDir = opt(refdir);
	const string &OutputFileName = opt(output);
	Dirize(TestDir1);
	Dirize(TestDir2);
	Dirize(RefDir);
	FILE *fOut = CreateStdioFile(OutputFileName);

	vector<string> Sets;
	ReadStringsFromFile(SetsFileName, Sets);

	const uint N = SIZE(Sets);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Running");
		const string &Set = Sets[i];
		string TestFileName1 = TestDir1 + Set;
		string TestFileName2 = TestDir2 + Set;
		string RefFileName = RefDir + Set;
		QScorer3 Q3;
		Q3.Run(TestFileName1, TestFileName2, RefFileName);
		Q3.ToTSV(fOut);
		}

	CloseStdioFile(fOut);
	}
