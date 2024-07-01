#include "muscle.h"
#include "qscorer.h"

void QScorer::Clear()
	{
	m_Test = 0;
	m_Ref = 0;
	m_RefAlignedColCount = 0;

	m_Labels.clear();
	m_RefSeqIndexes.clear();
	m_TestSeqIndexes.clear();
	m_RefSeqIndexToTestSeqIndex.clear();
	m_RefCols.clear();
	m_RefUngappedCounts.clear();

	m_PosToTestColVec.clear();
	m_PosToRefColVec.clear();

	m_TestColToPosVec.clear();
	m_RefColToPosVec.clear();

	m_RefColToTestColVec.clear();
	m_TestColToBestRefCol.clear();
	m_MaxFracts.clear();
	m_BestTestCols.clear();

	m_TotalPairs = 0;
	m_TotalCols = 0;

	m_CorrectPairs = 0;
	m_CorrectCols = 0;

	m_Q = 0;
	m_TC = 0;

	m_RefLabels.clear();
	m_RefLabelToSeqIndex.clear();
	m_TestColToCount.clear();

        m_RefSeqToSeqIndex.clear();
	}

void QScorer::InitRefLabels()
	{
	m_RefLabels.clear();
	m_RefLabelToSeqIndex.clear();

	const uint RefSeqCount = GetRefSeqCount();
	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const string Label = (const string) m_Ref->GetSeqName(RefSeqIndex);
		if (m_RefLabelToSeqIndex.find(Label) != m_RefLabelToSeqIndex.end())
			Die("Dupe ref label >%s", Label.c_str());

		m_RefLabels.push_back(Label);
		m_RefLabelToSeqIndex[Label] = RefSeqIndex;
		
		string Seq;
		m_Ref->GetUngappedSeqStr(RefSeqIndex, Seq);
		for (auto & c: Seq) c = toupper(c);
                if (m_RefSeqToSeqIndex.find(Seq) != m_RefSeqToSeqIndex.end())
		{
			uint RIndex = m_RefSeqToSeqIndex[Seq];
                        Die("Dupe ref seqs >%s %s", Label.c_str(),m_RefLabels[RIndex].c_str());
		}

                m_RefSeqToSeqIndex[Seq] = RefSeqIndex;



		}
	}
map<string, uint>::const_iterator QScorer::_FullScanRefSeqs(const string & TestSeq)
{
	map<string, uint>::const_iterator best = m_RefSeqToSeqIndex.end();
	uint min_dist =  UINT_MAX;
	for (map<string, uint>::const_iterator p = m_RefSeqToSeqIndex.begin(); p !=  m_RefSeqToSeqIndex.end(); p++)
	{
		uint dist = 0;
		const string & RefSeq = p->first;
		if (RefSeq.size() != TestSeq.size())
			continue;
		for (uint i = 0; i < RefSeq.size(); ++i)
			if (RefSeq[i] != TestSeq[i]) dist++;
		if (dist < min_dist)
		{
			min_dist = dist;
			best = p;
		}
	}
	if (min_dist < 4)
	{
		return best;
	}
	return m_RefSeqToSeqIndex.end();
}

void QScorer::InitRefToTest_BySequence()
	{
        const uint RefSeqCount = GetRefSeqCount();
        const uint TestSeqCount = GetTestSeqCount();

        m_TestSeqIndexToRefSeqIndex.clear();
        m_TestSeqIndexToRefSeqIndex.resize(TestSeqCount, UINT_MAX);

        m_RefSeqIndexToTestSeqIndex.clear();
        m_RefSeqIndexToTestSeqIndex.resize(RefSeqCount, UINT_MAX);
        for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
                {
                const string Label = (const string) m_Test->GetSeqName(TestSeqIndex);

                string Seq;
                m_Test->GetUngappedSeqStr(TestSeqIndex, Seq);
		for (auto & c: Seq) c = toupper(c);

                map<string, uint>::const_iterator p = m_RefSeqToSeqIndex.find(Seq);

                if (p == m_RefSeqToSeqIndex.end())
                {
			p = _FullScanRefSeqs(Seq);
			if (p == m_RefSeqToSeqIndex.end())
			{
Log("Seq mismatch >%s\n", Label.c_str());
                        continue;
			}
                }

                uint RefSeqIndex = p->second;
                asserta(RefSeqIndex < RefSeqCount);
                m_TestSeqIndexToRefSeqIndex[TestSeqIndex] = RefSeqIndex;
                if (m_RefSeqIndexToTestSeqIndex[RefSeqIndex] != UINT_MAX)
                        Warning("Ref label found twice in test MSA %s >%s",
                          m_Name.c_str(), Label.c_str());
                m_RefSeqIndexToTestSeqIndex[RefSeqIndex] = TestSeqIndex;

                m_Labels.push_back(Label);
                m_RefSeqIndexes.push_back(RefSeqIndex);
                m_TestSeqIndexes.push_back(TestSeqIndex);
                }
	}
void QScorer::InitRefToTest()
	{
    	if (opt_bysequence)
         {
         	InitRefToTest_BySequence();
         	return;
         }
	const uint RefSeqCount = GetRefSeqCount();
	const uint TestSeqCount = GetTestSeqCount();

	m_TestSeqIndexToRefSeqIndex.clear();
	m_TestSeqIndexToRefSeqIndex.resize(TestSeqCount, UINT_MAX);

	m_RefSeqIndexToTestSeqIndex.clear();
	m_RefSeqIndexToTestSeqIndex.resize(RefSeqCount, UINT_MAX);
	for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
		{
		const string Label = (const string) m_Test->GetSeqName(TestSeqIndex);
		map<string, uint>::const_iterator p = m_RefLabelToSeqIndex.find(Label);
		if (p == m_RefLabelToSeqIndex.end())
		{
Log("Label mismatch >%s\n", Label.c_str());
			continue;
		}
		uint RefSeqIndex = p->second;
		asserta(RefSeqIndex < RefSeqCount);
		m_TestSeqIndexToRefSeqIndex[TestSeqIndex] = RefSeqIndex;
		if (m_RefSeqIndexToTestSeqIndex[RefSeqIndex] != UINT_MAX)
			Warning("Ref label found twice in test MSA %s >%s",
			  m_Name.c_str(), Label.c_str());
		m_RefSeqIndexToTestSeqIndex[RefSeqIndex] = TestSeqIndex;

		m_Labels.push_back(Label);
		m_RefSeqIndexes.push_back(RefSeqIndex);
		m_TestSeqIndexes.push_back(TestSeqIndex);
		}
	}

void QScorer::InitColPosVecs1(uint i)
	{
	const uint TestColCount = GetTestColCount();
	const uint RefColCount = GetRefColCount();

	uint RefSeqIndex = m_RefSeqIndexes[i];
	uint TestSeqIndex = m_TestSeqIndexes[i];
	const string &Label = m_Labels[i];
	const string &TestLabel = m_Test->GetSeqName(TestSeqIndex);
	const string &RefLabel = m_Ref->GetSeqName(RefSeqIndex);
//	asserta(TestLabel == RefLabel);

	m_Ref->GetPosToCol(RefSeqIndex, m_PosToRefColVec[i]);
	m_Test->GetPosToCol(TestSeqIndex, m_PosToTestColVec[i]);
	const uint RefUngappedLength = SIZE(m_PosToRefColVec[i]);
	const uint TestUngappedLength = SIZE(m_PosToTestColVec[i]);
	if (RefUngappedLength != TestUngappedLength)
		Warning("%s >%s ref length %u, test length %u",
		  m_Name.c_str(), TestLabel.c_str(),
		  RefUngappedLength, TestUngappedLength);

	m_Ref->GetColToPos(RefSeqIndex, m_RefColToPosVec[i]);
	m_Test->GetColToPos(TestSeqIndex, m_TestColToPosVec[i]);

	asserta(SIZE(m_RefColToPosVec[i]) == RefColCount);
	asserta(SIZE(m_TestColToPosVec[i]) == TestColCount);

	const uint L = SIZE(m_PosToRefColVec[i]);
	const uint Lt = SIZE(m_PosToTestColVec[i]);
	if (L != Lt)
		Warning("Seq lengths differ ref=%u, test=%u >%s",
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
			char TestChar = m_Test->GetChar(TestSeqIndex, TestCol);
			char RefChar = m_Ref->GetChar(RefSeqIndex, RefCol);
			asserta(!isgap(TestChar) && !isgap(RefChar));
			char ut = toupper(TestChar);
			char ur = toupper(RefChar);
			if (ut != ur && ut != 'X' && ur != 'X')
				Warning("Sequences differ pos %u test %c ref %c >%s",
					Pos, TestChar, RefChar, Label.c_str());
			m_RefColToTestColVec[i][RefCol] = TestCol;
			}
		}
	}

void QScorer::InitColPosVecs()
	{
	const uint N = SIZE(m_TestSeqIndexes); // UINT_MAX if not in ref
	const uint NR = SIZE(m_RefLabels);
	asserta(N <= NR);
	if (N < NR && !m_MissingTestSeqOk)
		Die("%u missing sequences in test MSA %s", NR - N, m_Name.c_str());
	if (N == 0)
		{
		m_Q = -999;
		m_TC = -999;
		return;
		}

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
		InitColPosVecs1(i);
	}

void QScorer::InitRefCols()
	{
	const uint RefColCount = GetRefColCount();

	m_RefCols.clear();
	for (uint RefCol = 0; RefCol < RefColCount; ++RefCol)
		if (m_Ref->ColIsUpper(RefCol, m_MaxGapFract))
			m_RefCols.push_back(RefCol);
	}

void QScorer::InitRefUngappedCounts()
	{
	m_RefAlignedColCount = SIZE(m_RefCols);
	if (m_RefAlignedColCount == 0)
		Die("Qscorer: No upper case columns in ref %s", m_Name.c_str());

	m_RefUngappedCounts.clear();
	const uint N = SIZE(m_RefSeqIndexes);
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		{
		uint RefCol = m_RefCols[k];
		uint UngappedCount = 0;
		for (uint i = 0; i < N; ++i)
			{
			uint RefSeqIndex = m_RefSeqIndexes[i];
			char c = m_Ref->GetChar(RefSeqIndex, RefCol);
			if (!isgap(c))
				++UngappedCount;
			}
		m_RefUngappedCounts.push_back(UngappedCount);
		}
	}

void QScorer::DoRefCol(uint k)
	{
	uint RefCol = m_RefCols[k];
	const uint RefSeqCount = GetRefSeqCount();

	uint64 CorrectPairsCol = 0;
	vector<uint> TestColIndexesFound;
	const uint N = SIZE(m_RefSeqIndexes);
	uint TestLetterCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint TestCol = m_RefColToTestColVec[i][RefCol];
		if (TestCol != UINT_MAX)
			{
			++TestLetterCount;
			asserta(TestCol < SIZE(m_TestColToCount));
			if (m_TestColToCount[TestCol] == 0)
				TestColIndexesFound.push_back(TestCol);
			m_TestColToCount[TestCol] += 1;
			}
		}

	uint64 MaxCount = 0;
	uint BestTestCol = UINT_MAX;
	for (uint j = 0; j < SIZE(TestColIndexesFound); ++j)
		{
		uint TestCol = TestColIndexesFound[j];
		uint Count = m_TestColToCount[TestCol];
		asserta(Count > 0);
		uint Count64 = Count;
		if (Count > MaxCount)
			{
			MaxCount = Count;
			BestTestCol = TestCol;
			}
		CorrectPairsCol += (Count64*(Count64 - 1))/2;

	// Reset so all zero at start of next DoRefCol().
		m_TestColToCount[TestCol] = 0;
		}
	m_CorrectPairs += CorrectPairsCol;
	if (MaxCount <= TestLetterCount/2)
		BestTestCol = UINT_MAX;
	m_BestTestCols.push_back(BestTestCol);
	float MaxFract = float(MaxCount)/RefSeqCount;
	m_MaxFracts.push_back(MaxFract);

	asserta(k < SIZE(m_RefUngappedCounts));
	uint64 UngappedCount = m_RefUngappedCounts[k];
	uint64 UngappedPairCount = (UngappedCount*(UngappedCount - 1))/2;
	m_TotalPairs += UngappedPairCount;

	asserta(UngappedPairCount >= CorrectPairsCol);
	if (UngappedPairCount == CorrectPairsCol && CorrectPairsCol > 0)
		++m_CorrectCols;
	}

void QScorer::DoRefCols()
	{
	m_BestTestCols.clear();
	m_MaxFracts.clear();
	m_TestColToCount.clear();

	const uint TestColCount = GetTestColCount();
	m_TestColToCount.resize(TestColCount, 0);

	m_CorrectPairs = 0;
	m_CorrectCols = 0;
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		DoRefCol(k);
	}

void QScorer::SetTestColToBestRefCol()
	{
	const uint TestColCount = GetTestColCount();
	m_TestColToBestRefCol.resize(TestColCount, UINT_MAX);
	for (uint k = 0; k < m_RefAlignedColCount; ++k)
		{
		uint RefCol = m_RefCols[k];
		uint BestTestCol = m_BestTestCols[k];
		if (RefCol == UINT_MAX || BestTestCol == UINT_MAX)
			continue;

		asserta(BestTestCol < SIZE(m_TestColToBestRefCol));
		m_TestColToBestRefCol[BestTestCol] = RefCol;
		}
	}

void QScorer::Run(const string &Name, const MultiSequence &Test, const MultiSequence &Ref)
	{
	MSA &msaTest = *new MSA;
	MSA &msaRef = *new MSA;

	msaTest.FromMultiSequence(Test);
	msaRef.FromMultiSequence(Ref);
	Run(Name, msaTest, msaRef);

	//delete &msaTest; // crashes, just leak it
	//delete &msaRef;
	}

void QScorer::Run(const string &Name, const MSA &Test, const MSA &Ref)
	{
	Clear();

	m_Name = Name;
	m_Q = 0;
	m_TC = 0;
	m_Test = &Test;
	m_Ref = &Ref;

	const uint TestSeqCount = Test.GetSeqCount();
	const uint RefSeqCount = Ref.GetSeqCount();

	const uint TestColCount = Test.GetColCount();
	const uint RefColCount = Ref.GetColCount();

	InitRefLabels();
	asserta(SIZE(m_RefLabels) > 0);
	InitRefToTest();
	InitColPosVecs();
	if (SIZE(m_RefSeqIndexes) <= 1)
		{
		m_Q = -1.0;
		m_TC = -1.0;
		return;
		}
	InitRefCols();
	InitRefUngappedCounts();
	DoRefCols();
	SetTestColToBestRefCol();

	if (m_TotalPairs == 0)
		Die("m_TotalPairs=0 (%s)", m_Name.c_str());
	asserta(m_RefAlignedColCount > 0);
	m_Q = float(m_CorrectPairs)/float(m_TotalPairs);
	m_TC = float(m_CorrectCols)/float(m_RefAlignedColCount);
	}

void QScorer::UpdateRefLetterCountsCol(uint k,
  vector<vector<uint> > &LetterCountsVec) const
	{
	asserta(k < SIZE(m_RefCols));
	asserta(k < SIZE(m_BestTestCols));

	uint RefCol = m_RefCols[k];
	uint BestTestCol = m_BestTestCols[k];

	const uint N = SIZE(m_RefSeqIndexes);
	asserta(SIZE(m_TestSeqIndexes) == N);
	for (uint i = 0; i < N; ++i)
		{
		uint Pos = m_RefColToPosVec[i][RefCol];
		if (Pos == UINT_MAX)
			continue;
		uint TestCol = m_PosToTestColVec[i][Pos];
		if (TestCol == BestTestCol)
			{
			uint RefSeqIndex = m_RefSeqIndexes[i];
			LetterCountsVec[RefSeqIndex][RefCol] += 1;
			}
		}
	}

/***
LetterCountsVec[RefSeqIndex][RefColIndex] is the number of times
this position appears in the best-match test column, will be zero
or one on first pass, can be incremented by calling again with
more test MSAs with the same ref MSA.
***/
void QScorer::UpdateRefLetterCounts(vector<vector<uint> > &LetterCountsVec) const
	{
	const uint RefSeqCount = GetRefSeqCount();
	const uint RefColCount = GetRefColCount();

	if (LetterCountsVec.empty())
		{
		LetterCountsVec.clear();
		LetterCountsVec.resize(RefSeqCount);
		for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
			LetterCountsVec[RefSeqIndex].resize(RefColCount, 0);
		}
	else
		{
		asserta(RefSeqCount > 0);
		asserta(SIZE(LetterCountsVec) == RefSeqCount);
		asserta(SIZE(LetterCountsVec[0]) == RefColCount);
		}

	const uint K = SIZE(m_RefCols);
	for (uint k = 0; k < K; ++k)
		UpdateRefLetterCountsCol(k, LetterCountsVec);
	}
