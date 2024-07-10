#include "muscle.h"
#include "qscorer.h"

void AlignMSAsByCols(const MSA &X, const MSA &Y,
  const vector<uint> &aColsX, const vector<uint> &aColsY,
  string &Path, MSA &X2, MSA &Y2);

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

	MSA TestSubset;
	MSA RefSubset;
	MSAFromSeqSubset(Test, TestSeqIndexes.data(), SIZE(TestSeqIndexes), TestSubset);
	MSAFromSeqSubset(Ref, RefSeqIndexes.data(), SIZE(RefSeqIndexes), RefSubset);

	MSA Test2;
	MSA Ref2;
	string Path;
	AlignMSAsByCols(TestSubset, RefSubset, TestCols, RefCols, Path, Test2, Ref2);
	const uint ColCount2 = Test2.GetColCount();
	asserta(Ref2.GetColCount() == ColCount2);

	const uint PL = SIZE(Path);
	uint FirstM = UINT_MAX;
	uint LastM = UINT_MAX;
	for (uint i = 0; i < PL; ++i)
		{
		if (Path[i] == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = i;
			LastM = i;
			}
		}

	const uint TestSeqCount = Test2.GetSeqCount();
	for (uint i = 0; i < TestSeqCount; ++i)
		{
		const char *Label = Test2.GetLabel(i);
		string Row;
		Test2.GetRowStr(i, Row);
		string RowM;
		for (uint j = FirstM; j <= LastM; ++j)
			RowM += Row[j];
		Log("%s  >%s\n", RowM.c_str(), Label);
		}
	Log("\n");
	const uint RefSeqCount = Ref2.GetSeqCount();
	for (uint i = 0; i < RefSeqCount; ++i)
		{
		const char *Label = Ref2.GetLabel(i);
		string Row;
		Ref2.GetRowStr(i, Row);
		string RowM;
		for (uint j = FirstM; j <= LastM; ++j)
			RowM += Row[j];
		Log("%s  >%s\n", RowM.c_str(), Label);
		}
	}
