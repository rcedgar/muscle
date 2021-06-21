#include "muscle.h"
#include "qscorer.h"
#include "qscorer3.h"

void QScorer3::TransQPair(uint Indexi, uint Indexj)
	{
	const string &Labeli = m_QS1.m_Labels[Indexi];
	const string &Labelj = m_QS1.m_Labels[Indexj];

	const uint RefSeqIndexi = m_QS1.m_RefSeqIndexes[Indexi];
	const uint RefSeqIndexj = m_QS1.m_RefSeqIndexes[Indexj];
	asserta(string(m_Ref.GetSeqName(RefSeqIndexi)) == Labeli);
	asserta(string(m_Ref.GetSeqName(RefSeqIndexj)) == Labelj);

	uint Indexi2 = m_Indexes2[Indexi];
	uint Indexj2 = m_Indexes2[Indexj];
	asserta(m_QS2.m_Labels[Indexi2] == Labeli);
	asserta(m_QS2.m_Labels[Indexj2] == Labelj);
	asserta(m_QS2.m_RefSeqIndexes[Indexi2] == RefSeqIndexi);
	asserta(m_QS2.m_RefSeqIndexes[Indexj2] == RefSeqIndexj);

	const vector<uint> &TestColToPosi1 = m_QS1.m_TestColToPosVec[Indexi];
	const vector<uint> &TestColToPosi2 = m_QS2.m_TestColToPosVec[Indexi2];

	const vector<uint> &TestColToPosj1 = m_QS1.m_TestColToPosVec[Indexj];
	const vector<uint> &TestColToPosj2 = m_QS2.m_TestColToPosVec[Indexj2];

	const vector<uint> &RefColToPosi1 = m_QS1.m_RefColToPosVec[Indexi];
	const vector<uint> &RefColToPosi2 = m_QS2.m_RefColToPosVec[Indexi2];

	const vector<uint> &RefColToPosj1 = m_QS1.m_RefColToPosVec[Indexj];
	const vector<uint> &RefColToPosj2 = m_QS2.m_RefColToPosVec[Indexj2];

	const vector<uint> &PosToTestColi1 = m_QS1.m_PosToTestColVec[Indexi];
	const vector<uint> &PosToTestColi2 = m_QS2.m_PosToTestColVec[Indexi2];

	const vector<uint> &PosToTestColj1 = m_QS1.m_PosToTestColVec[Indexj];
	const vector<uint> &PosToTestColj2 = m_QS2.m_PosToTestColVec[Indexj2];

	const uint RefColCount = SIZE(RefColToPosi1);
	asserta(SIZE(RefColToPosi2) == RefColCount);
	asserta(SIZE(RefColToPosj1) == RefColCount);
	asserta(SIZE(RefColToPosj2) == RefColCount);

	const uint Li = SIZE(PosToTestColi1);
	const uint Lj = SIZE(PosToTestColj1);

	asserta(SIZE(PosToTestColi2) == Li);
	asserta(SIZE(PosToTestColj2) == Lj);

	const vector<uint> &RefCols = *m_RefCols;
	const uint RefAlignedColCount = m_QS1.m_RefAlignedColCount;
	asserta(m_QS2.m_RefAlignedColCount == RefAlignedColCount);

	uint CorrectColCount1 = 0;
	uint CorrectColCount2 = 0;
	vector<uint> Posis;
	vector<uint> Posjs;
	for (uint k = 0; k < RefAlignedColCount; ++k)
		{
		uint RefCol = RefCols[k];

		uint Posi = RefColToPosi1[RefCol];
		uint Posi2 = RefColToPosi2[RefCol];
		if (Posi == UINT_MAX || Posi2 == UINT_MAX)
			continue;
		asserta(Posi2 == Posi);

		uint TestColi1 = PosToTestColi1[Posi];
		uint TestColi2 = PosToTestColi2[Posi];

		uint Posj = RefColToPosj1[RefCol];
		uint Posj2 = m_QS2.m_RefColToPosVec[Indexj2][RefCol];
		if (Posj == UINT_MAX || Posj2 == UINT_MAX)
			continue;
		asserta(Posj2 == Posj);
		Posis.push_back(Posi);
		Posjs.push_back(Posj);

		uint TestColj1 = PosToTestColj1[Posj];
		uint TestColj2 = PosToTestColj2[Posj];
				
		if (TestColi1 == TestColj1)
			++CorrectColCount1;
		if (TestColi2 == TestColj2)
			++CorrectColCount2;
		}

	const uint TestAlignedPosCount = SIZE(Posis);
	asserta(SIZE(Posjs) == TestAlignedPosCount);

	uint SameColCount = 0;
	for (uint k = 0; k < TestAlignedPosCount; ++k)
		{
		uint Posi = Posis[k];

		uint TestCol1 = PosToTestColi1[Posi];
		uint Posj1 = TestColToPosj1[TestCol1];

		uint TestCol2 = PosToTestColi2[Posi];
		uint Posj2 = TestColToPosj2[TestCol2];

		if (Posj1 == Posj2)
			++SameColCount;
		}

	for (uint k = 0; k < TestAlignedPosCount; ++k)
		{
		uint Posj = Posjs[k];

		uint TestCol1 = PosToTestColj1[Posj];
		uint Posi1 = TestColToPosi1[TestCol1];

		uint TestCol2 = PosToTestColj2[Posj];
		uint Posi2 = TestColToPosi2[TestCol2];

		if (Posi1 == Posi2)
			++SameColCount;
		}

	float Q1 = float(CorrectColCount1)/RefAlignedColCount;
	float Q2 = float(CorrectColCount2)/RefAlignedColCount;;
	float PWC = float(SameColCount)/(2*TestAlignedPosCount);

	m_Pairs.push_back(pair<uint, uint>(Indexi, Indexj));
	m_PairIndexToQ1.push_back(Q1);
	m_PairIndexToQ2.push_back(Q2);
	m_PairIndexToPWC.push_back(PWC);
	}

void QScorer3::TransQ()
	{
	const vector<string> &Labels = m_QS1.m_Labels;
	const uint N = SIZE(Labels);

	asserta(m_RefCols != 0);
	const vector<uint> &RefCols = *m_RefCols;

	m_Pairs.clear();
	m_PairIndexToQ1.clear();
	m_PairIndexToQ2.clear();
	m_PairIndexToPWC.clear();
	for (uint Indexi = 0; Indexi < N; ++Indexi)
		for (uint Indexj = Indexi + 1; Indexj < N; ++Indexj)
			TransQPair(Indexi, Indexj);
	}
