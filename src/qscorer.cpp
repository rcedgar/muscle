#include "muscle.h"
#include "qscorer.h"
#include "quarts.h"

void QScorer::Run(const string &TestName, const string &RefName,
  const MSA &Test, const MSA &Ref)
	{
	m_TestName = TestName;
	m_RefName = RefName;

	m_Test = &Test;
	m_Ref = &Ref;

	const uint TestSeqCount = Test.GetSeqCount();
	const uint RefSeqCount = Ref.GetSeqCount();

	const uint TestColCount = Test.GetColCount();
	const uint RefColCount = Ref.GetColCount();

	vector<string> RefLabels;
	map<string, uint> RefLabelToSeqIndex;
	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const string Label = (const string) Ref.GetSeqName(RefSeqIndex);
		if (RefLabelToSeqIndex.find(Label) != RefLabelToSeqIndex.end())
			Die("Dupe ref label >%s", Label.c_str());

		RefLabels.push_back(Label);
		RefLabelToSeqIndex[Label] = RefSeqIndex;
		}

	m_RefSeqIndexes.clear();
	m_TestSeqIndexes.clear();
	m_Labels.clear();

	m_RefSeqIndexToTestSeqIndex.clear();
	m_RefSeqIndexToTestSeqIndex.resize(RefSeqCount, UINT_MAX);
	for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
		{
		const string Label = (const string) Test.GetSeqName(TestSeqIndex);
		map<string, uint>::const_iterator p = RefLabelToSeqIndex.find(Label);
		if (p == RefLabelToSeqIndex.end())
			continue;
		uint RefSeqIndex = p->second;
		asserta(RefSeqIndex < RefSeqCount);
		if (m_RefSeqIndexToTestSeqIndex[RefSeqIndex] != UINT_MAX)
			Die("Ref label found twice in test MSA >%s", Label.c_str());
		m_RefSeqIndexToTestSeqIndex[RefSeqIndex] = TestSeqIndex;

		m_Labels.push_back(Label);
		m_RefSeqIndexes.push_back(RefSeqIndex);
		m_TestSeqIndexes.push_back(TestSeqIndex);
		}

	const uint N = SIZE(m_RefSeqIndexes);
	if (N == 0)
		Die("No ref labels found in test MSA");

	m_PosToTestColVec.clear();
	m_PosToRefColVec.clear();

	m_TestColToPosVec.clear();
	m_RefColToPosVec.clear();

	m_RefColToTestColVec.clear();

	m_PosToTestColVec.resize(N);
	m_PosToRefColVec.resize(N);

	m_TestColToPosVec.resize(N);
	m_RefColToPosVec.resize(N);

	m_RefColToTestColVec.resize(N);

	for (uint i = 0; i < N; ++i)
		{
		uint RefSeqIndex = m_RefSeqIndexes[i];
		uint TestSeqIndex = m_TestSeqIndexes[i];
		const string &Label = m_Labels[i];
		const string &TestLabel = Test.GetSeqName(TestSeqIndex);
		const string &RefLabel = Ref.GetSeqName(RefSeqIndex);
		asserta(TestLabel == RefLabel);

		Ref.GetPosToCol(RefSeqIndex, m_PosToRefColVec[i]);
		Test.GetPosToCol(TestSeqIndex, m_PosToTestColVec[i]);
		const uint RefUngappedLength = SIZE(m_PosToRefColVec[i]);
		const uint TestUngappedLength = SIZE(m_PosToTestColVec[i]);
		asserta(RefUngappedLength == TestUngappedLength);

		Ref.GetColToPos(RefSeqIndex, m_RefColToPosVec[i]);
		Test.GetColToPos(TestSeqIndex, m_TestColToPosVec[i]);

		asserta(SIZE(m_RefColToPosVec[i]) == RefColCount);
		asserta(SIZE(m_TestColToPosVec[i]) == TestColCount);

		const uint L = SIZE(m_PosToRefColVec[i]);
		const uint Lt = SIZE(m_PosToTestColVec[i]);
		if (L != Lt)
			Die("Seq lengths differ ref=%u, test=%u >%s",
			  L, Lt, Label.c_str());

		m_RefColToTestColVec[i].resize(RefColCount, UINT_MAX);
		for (uint RefCol = 0; RefCol < RefColCount; ++RefCol)
			{
			uint Pos = m_RefColToPosVec[i][RefCol];
			if (Pos == UINT_MAX)
				m_RefColToTestColVec[i][RefCol] = UINT_MAX;
			else
				{
				asserta(Pos < SIZE(m_PosToTestColVec[i]));
				uint TestCol = m_PosToTestColVec[i][Pos];
				asserta(TestCol < TestColCount);
				char TestChar = Test.GetChar(TestSeqIndex, TestCol);
				char RefChar = Ref.GetChar(RefSeqIndex, RefCol);
				asserta(!isgap(TestChar) && !isgap(RefChar));
				if (toupper(TestChar) != toupper(RefChar))
					Die("Sequences differ pos %u test %c ref %c >%s",
					  Pos, TestChar, RefChar, Label.c_str());
				m_RefColToTestColVec[i][RefCol] = TestCol;
				}
			}
		}

	m_RefCols.clear();
	for (uint RefCol = 0; RefCol < RefColCount; ++RefCol)
		if (Ref.ColIsUpper(RefCol))
			m_RefCols.push_back(RefCol);

	m_RefAlignedColCount = SIZE(m_RefCols);
	if (m_RefAlignedColCount == 0)
		Die("No upper case columns in ref %s", RefName.c_str());

	m_RefUngappedCounts.clear();
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		{
		uint RefCol = m_RefCols[k];
		uint UngappedCount = 0;
		for (uint i = 0; i < N; ++i)
			{
			uint RefSeqIndex = m_RefSeqIndexes[i];
			char c = Ref.GetChar(RefSeqIndex, RefCol);
			if (!isgap(c))
				++UngappedCount;
			}
		m_RefUngappedCounts.push_back(UngappedCount);
		}

	vector<uint> TestColToCount(TestColCount, 0);
	m_BestTestCols.clear();
	m_MaxFracts.clear();

	m_CorrectPairs = 0;
	m_CorrectCols = 0;
	m_CorrectCols90 = 0;
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		{
		uint RefCol = m_RefCols[k];

		uint64 CorrectPairsCol = 0;
		vector<uint> TestColIndexesFound;
		for (uint i = 0; i < N; ++i)
			{
			uint TestCol = m_RefColToTestColVec[i][RefCol];
			if (TestCol != UINT_MAX)
				{
				asserta(TestCol < SIZE(TestColToCount));
				if (TestColToCount[TestCol] == 0)
					TestColIndexesFound.push_back(TestCol);
				TestColToCount[TestCol] += 1;
				}
			}

		uint64 MaxCount = 0;
		uint BestTestCol = UINT_MAX;
		for (uint j = 0; j < SIZE(TestColIndexesFound); ++j)
			{
			uint TestCol = TestColIndexesFound[j];
			uint Count = TestColToCount[TestCol];
			asserta(Count > 0);
			TestColToCount[TestCol] = 0;
			uint Count64 = Count;
			if (Count > MaxCount)
				{
				MaxCount = Count;
				BestTestCol = TestCol;
				}
			CorrectPairsCol += (Count64*(Count64 - 1))/2;
			}
		m_CorrectPairs += CorrectPairsCol;
		m_BestTestCols.push_back(BestTestCol);
		float MaxFract = float(MaxCount)/RefSeqCount;
		m_MaxFracts.push_back(MaxFract);

		asserta(k < SIZE(m_RefUngappedCounts));
		uint64 UngappedCount = m_RefUngappedCounts[k];
		uint64 UngappedPairCount = (UngappedCount*(UngappedCount - 1))/2;
		m_TotalPairs += UngappedPairCount;

		asserta(UngappedPairCount >= CorrectPairsCol);
		if (UngappedPairCount == CorrectPairsCol)
			++m_CorrectCols;
		if (UngappedPairCount >= uint((CorrectPairsCol*90.0 + 0.5)/100.0))
			++m_CorrectCols90;
		}

	m_TestColToBestRefCol.resize(TestColCount, UINT_MAX);
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		{
		uint RefCol = m_RefCols[k];
		uint BestTestCol = m_BestTestCols[k];
		if (RefCol != UINT_MAX)
			{
			asserta(BestTestCol < SIZE(m_TestColToBestRefCol));
			m_TestColToBestRefCol[BestTestCol] = RefCol;
			}
		}

	QuartsFloat QF;
	GetQuartsFloat(m_MaxFracts, QF);

	m_Q = float(m_CorrectPairs)/float(m_TotalPairs);
	m_TC = float(m_CorrectCols)/float(m_RefAlignedColCount);
	m_TC90 = float(m_CorrectCols90)/float(m_RefAlignedColCount);
	m_CF = QF.Avg;
	m_Acc = (m_Q + m_TC90)/2;
	}

void QScorer::Report1(FILE *f, uint Index) const
	{
	if (f == 0)
		return;

	asserta(m_Test != 0);
	const uint TestColCount = m_Test->GetColCount();

	asserta(Index < SIZE(m_TestSeqIndexes));
	asserta(Index < SIZE(m_RefSeqIndexes));
	const uint RefSeqIndex = m_RefSeqIndexes[Index];
	const uint TestSeqIndex = m_TestSeqIndexes[Index];

	const string &TestLabel = m_Test->GetSeqName(TestSeqIndex);
	const string &RefLabel = m_Ref->GetSeqName(RefSeqIndex);
	asserta(TestLabel == RefLabel);

	asserta(SIZE(m_TestColToBestRefCol) == TestColCount);

	const vector<uint> &PosToRefCol = m_PosToRefColVec[Index];
	const vector<uint> &PosToTestCol = m_PosToTestColVec[Index];
	const uint RefUngappedLength = SIZE(PosToRefCol);
	const uint UngappedLength = SIZE(PosToTestCol);
	if (RefUngappedLength != UngappedLength)
		Die("i %u test %u ref %u RefL %u TestL %u",
		  Index, TestSeqIndex, RefSeqIndex, RefUngappedLength, UngappedLength);

	const vector<uint> &TestColToPos = m_TestColToPosVec[TestSeqIndex];
	asserta(SIZE(TestColToPos) == TestColCount);

	string Row;
	for (uint TestCol = 0; TestCol < TestColCount; ++TestCol)
		{
		uint BestRefCol = m_TestColToBestRefCol[TestCol];
		if (BestRefCol == UINT_MAX)
			continue;

		char TestChar = m_Test->GetChar(TestSeqIndex, TestCol);
		if (isgap(TestChar))
			{
			Row += '-';
			continue;
			}

		asserta(TestCol < SIZE(TestColToPos));
		uint Pos = TestColToPos[TestCol];
		asserta(Pos < UngappedLength);

		asserta(Pos < SIZE(PosToRefCol));
		uint RefCol = PosToRefCol[Pos];
		asserta(RefCol != UINT_MAX);

		if (RefCol == BestRefCol)
			Row += toupper(TestChar);
		else
			Row += tolower(TestChar);
		}
	fprintf(f, "%s\n", Row.c_str());
	}

void QScorer::Report(FILE *f) const
	{
	if (f == 0)
		return;

	const uint N = SIZE(m_RefSeqIndexes);
	asserta(N > 0);
	for (uint i = 0; i < N; ++i)
		Report1(f, i);
	}
