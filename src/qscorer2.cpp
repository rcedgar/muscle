#include "muscle.h"
#include "qscorer2.h"

void QScorer2::StripGaps(const string &Row, string &Seq) const
	{
	Seq.resize(0);
	for (string::const_iterator iter = Row.begin();
	  iter != Row.end(); ++iter)
		{
		char c = *iter;
		if (!isgap(c))
			Seq += c;
		}
	}

void QScorer2::GetPosToCol(const string &Row,
  vector<uint> &PosToCol) const
	{
	PosToCol.resize(0);
	uint ColCount = SIZE(Row);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Row[Col];
		if (!isgap(c))
			PosToCol.push_back(Col);
		}
	}

double QScorer2::GetQ(const string &T1, const string &T2,
  const string &R1, const string &R2) const
	{
	const uint NT = SIZE(T1);
	asserta(SIZE(T2) == NT);

	const uint NR = SIZE(R1);
	asserta(SIZE(R2) == NR);

#if DEBUG
	{
	string UngappedT1;
	string UngappedT2;
	string UngappedR1;
	string UngappedR2;
	StripGaps(T1, UngappedT1);
	StripGaps(T2, UngappedT2);
	StripGaps(R1, UngappedR1);
	StripGaps(R2, UngappedR2);
	asserta(UngappedT1 == UngappedR1);
	asserta(UngappedT2 == UngappedR2);
	}
#endif

	vector<uint> PosToColT1;
	vector<uint> PosToColT2;
	GetPosToCol(T1, PosToColT1);
	GetPosToCol(T2, PosToColT2);

	uint RefPos1 = 0;
	uint RefPos2 = 0;
	uint RefAlignedCount = 0;
	uint CorrectCount = 0;
	for (uint RefCol = 0; RefCol < NR; ++RefCol)
		{
		char rc1 = R1[RefCol];
		char rc2 = R2[RefCol];
		bool gap1 = isgap(rc1);
		bool gap2 = isgap(rc2);
		if (!gap1 && !gap2)
			{
			asserta(RefPos1 < SIZE(PosToColT1));
			asserta(RefPos2 < SIZE(PosToColT2));
			uint TestCol1 = PosToColT1[RefPos1];
			uint TestCol2 = PosToColT2[RefPos2];
			char tc1 = T1[TestCol1];
			char tc2 = T2[TestCol2];
			asserta(toupper(rc1) == toupper(tc1));
			asserta(toupper(rc2) == toupper(tc2));
			++RefAlignedCount;
			if (TestCol1 == TestCol2)
				++CorrectCount;
			}
		if (!gap1)
			++RefPos1;
		if (!gap2)
			++RefPos2;
		}
	double Q = double(CorrectCount)/RefAlignedCount;
	return Q;
	}

double QScorer2::Run(const MultiSequence &Test, const MultiSequence &Ref)
	{
	MSA *msaTest = new MSA;
	msaTest->FromMultiSequence(Test);
	double Q = Run(*msaTest, Ref);
	//delete msaTest; // crashes, just leak it
	return Q;
	}

double QScorer2::Run(const MSA &Test, const MultiSequence &Ref)
	{
	m_Test = &Test;
	m_Ref = &Ref;

	vector<string> Labels;
	map<string, uint> TestLabelToSeqIndex;
	Test.GetLabelToSeqIndex(Labels, TestLabelToSeqIndex);

	const uint RefSeqCount = Ref.GetSeqCount();
	asserta(RefSeqCount%2 == 0);
	const uint RefPairCount = RefSeqCount/2;
	double SumQ = 0;
	for (uint RefPairIndex = 0; RefPairIndex < RefPairCount;
	  ++RefPairIndex)
		{
		uint RefSeqIndex1 = RefPairIndex*2;
		uint RefSeqIndex2 = RefSeqIndex1 + 1;
		const string &Label1 = Ref.GetLabel(RefSeqIndex1);
		const string &Label2 = Ref.GetLabel(RefSeqIndex2);
		uint ColCount = Ref.GetSeqLength(RefSeqIndex1);
		asserta(Ref.GetSeqLength(RefSeqIndex2) == ColCount);

		map<string, uint>::const_iterator iter1 =
		  TestLabelToSeqIndex.find(Label1);
		map<string, uint>::const_iterator iter2 =
		  TestLabelToSeqIndex.find(Label2);
		asserta(iter1 != TestLabelToSeqIndex.end());
		asserta(iter2 != TestLabelToSeqIndex.end());

		const uint TestSeqIndex1 = iter1->second;
		const uint TestSeqIndex2 = iter2->second;

		string TestRow1;
		string TestRow2;
		Test.GetRowStr(TestSeqIndex1, TestRow1);
		Test.GetRowStr(TestSeqIndex2, TestRow2);

		string RefRow1;
		string RefRow2;
		Ref.GetSeqStr(RefSeqIndex1, RefRow1);
		Ref.GetSeqStr(RefSeqIndex2, RefRow2);

		double Q = GetQ(TestRow1, TestRow2, RefRow1, RefRow2);
		SumQ += Q;
		}
	double AvgQ = SumQ/RefPairCount;
	return AvgQ;
	}
