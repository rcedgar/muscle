#include "muscle.h"
#include "qscorer.h"

void AlignMSAsByCols(const MSA &X, const MSA &Y,
  const vector<uint> &ColsX, const vector<uint> &ColsY,
  string &Path, vector<uint> &MergeMap, MSA &X2, MSA &Y2);

static char GetQChar(double Q)
	{
	asserta(Q >= 0 && Q <= 1);
	if (Q == 1)
		return '|';
	else if (Q >= 0.9)
		return ':';
	else if (Q >= 0.5)
		return '.';
	else if (Q > 0)
		return '@';
	return '*';
	}

void cmd_cmp_ref_msas()
	{
	const string TestFileName = g_Arg1;
	const string RefFileName = opt(ref);

	double MaxGapFract = optd(max_gap_fract, 1.0);

	string Name;
	GetBaseName(TestFileName.c_str(), Name);

	MSA Test;
	MSA Ref;
	Test.FromFASTAFile_PreserveCase(TestFileName);
	Ref.FromFASTAFile_PreserveCase(RefFileName);

	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;
	QS.CmpRefMSAs(Name, Test, Ref);

	const uint NC = QS.m_RefMSAs_ComparedColCount;
	ProgressLog("@CMP_REF_MSAs test=%s ref=%s name=%s cols=%u Q=%.4f\n",
	  TestFileName.c_str(), RefFileName.c_str(), Name.c_str(), NC, QS.m_RefMSAs_Q);
	if (NC == 0)
		return;
	const vector<double> &Qs = QS.m_RefMSAs_ColQs;
	const vector<uint> &TestCols = QS.m_RefMSAs_TestCols;
	const vector<uint> &RefCols = QS.m_RefMSAs_RefCols;
	asserta(SIZE(TestCols) == NC);
	asserta(SIZE(RefCols) == NC);
	asserta(SIZE(Qs) == NC);
	const vector<uint> &TestSeqIndexes = QS.m_TestSeqIndexes;
	const vector<uint> &RefSeqIndexes = QS.m_RefSeqIndexes;
	const uint NS = SIZE(TestSeqIndexes);
	asserta(SIZE(RefSeqIndexes) == NS);
	for (uint i = 0; i < NC; ++i)
		{
		uint TestCol = TestCols[i];
		uint RefCol = RefCols[i];
		string TestColStr;
		string RefColStr;
		for (uint j = 0; j < NS; ++j)
			{
			uint TestSeqIndex = TestSeqIndexes[j];
			uint RefSeqIndex = RefSeqIndexes[j];
			TestColStr += Test.GetChar(TestSeqIndex, TestCol);
			RefColStr += Ref.GetChar(RefSeqIndex, RefCol);
			}
		double Q = Qs[i];
		Log("%s  %s  %6.4f\n", TestColStr.c_str(), RefColStr.c_str(), Q);
		}

	Log("\n");
	const uint aNC = SIZE(TestCols);
	asserta(aNC > 0);
	vector<uint> FixedTestCols;
	vector<uint> FixedRefCols;
	vector<double> FixedQs;
	FixedTestCols.push_back(TestCols[0]);
	FixedRefCols.push_back(RefCols[0]);
	FixedQs.push_back(Qs[0]);
	for (uint i = 1; i < aNC; ++i)
		{
		uint ColX = TestCols[i];
		uint ColY = RefCols[i];
		if (ColX > FixedTestCols.back() && ColY > FixedRefCols.back())
			{
			FixedTestCols.push_back(ColX);
			FixedRefCols.push_back(ColY);
			FixedQs.push_back(Qs[i]);
			}
		}
	const uint NF = SIZE(FixedTestCols);
	asserta(SIZE(FixedRefCols) == NF);
	asserta(NF > 0);

	MSA TestSubset;
	MSA RefSubset;
	MSAFromSeqSubset(Test, TestSeqIndexes.data(), SIZE(TestSeqIndexes), TestSubset);
	MSAFromSeqSubset(Ref, RefSeqIndexes.data(), SIZE(RefSeqIndexes), RefSubset);

	MSA Test2;
	MSA Ref2;
	string Path;
	vector<uint> MergeMap;
	AlignMSAsByCols(TestSubset, RefSubset,
	  FixedTestCols, FixedRefCols, Path, MergeMap, Test2, Ref2);
	const uint ColCount2 = SIZE(Path);
	asserta(Test2.GetColCount() == ColCount2);
	asserta(Ref2.GetColCount() == ColCount2);
	asserta(SIZE(MergeMap) == NF);

	vector<bool> AllGaps;
	for (uint i = 0; i < ColCount2; ++i)
		AllGaps.push_back(Test2.IsGapColumn(i) && Ref2.IsGapColumn(i));

	uint FirstM = UINT_MAX;
	uint LastM = UINT_MAX;
	for (uint i = 0; i < ColCount2; ++i)
		{
		if (Path[i] == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = i;
			LastM = i;
			}
		}

	string Annot(ColCount2, '_');
	for (uint i = 0; i < NF; ++i)
		{
		uint k = MergeMap[i];
		double Q = FixedQs[i];
		asserta(k < SIZE(Annot));
		Annot[k] = GetQChar(Q);
		}

	const uint TestSeqCount = Test2.GetSeqCount();
	for (uint i = 0; i < TestSeqCount; ++i)
		{
		const char *Label = Test2.GetLabel(i);
		const char *RefLabel = Ref2.GetLabel(i);
		string Row;
		Test2.GetRowStr(i, Row);
		string RowM;
		for (uint j = FirstM; j <= LastM; ++j)
			if (!AllGaps[j])
				RowM += Row[j];
		Log("%s  >%s (=%s)\n", RowM.c_str(), Label, RefLabel);
		}
	Log("\n");
	string AnnotM;
	for (uint j = FirstM; j <= LastM; ++j)
		if (!AllGaps[j])
			AnnotM += Annot[j];
	Log("%s\n", AnnotM.c_str());
	Log("\n");
	const uint RefSeqCount = Ref2.GetSeqCount();
	for (uint i = 0; i < RefSeqCount; ++i)
		{
		const char *Label = Ref2.GetLabel(i);
		const char *TestLabel = Test2.GetLabel(i);
		string Row;
		Ref2.GetRowStr(i, Row);
		string RowM;
		for (uint j = FirstM; j <= LastM; ++j)
			if (!AllGaps[j])
				RowM += Row[j];
		Log("%s  >%s\n", RowM.c_str(), Label, TestLabel);
		}
	}
