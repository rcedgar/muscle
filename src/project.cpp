#include "muscle.h"

MultiSequence* MultiSequence::Project(const vector<uint> &SeqIndexes)
	{
	set<int> IndexSet;
	for (uint i = 0; i < SIZE(SeqIndexes); ++i)
		{
		int Index = (int) SeqIndexes[i];
		IndexSet.insert(Index);
		}
	return Project(IndexSet);
	}

// Create new MSA from subset of sequences, deleting
// all-blank columns
MultiSequence* MultiSequence::Project(const set<int>& SeqIndexes)
	{
	MultiSequence* SubsetMSA = new MultiSequence;
	vector<Sequence *> NewSeqs;
	uint OldColCount = UINT_MAX;
	for (set<int>::const_iterator p = SeqIndexes.begin();
	  p != SeqIndexes.end(); ++p)
		{
		int OldSeqIndex = *p;
		const Sequence *OldSeq = GetSequence(OldSeqIndex);
		const string &Label = OldSeq->m_Label;
		uint L = OldSeq->GetLength();
		if (OldColCount == UINT_MAX)
			OldColCount = L;
		else
			asserta(L == OldColCount);

		Sequence *NewSeq = NewSequence(); 
		NewSeq->m_Label = Label;
		NewSeq->m_CharVec.reserve(L);
		SubsetMSA->m_Seqs.push_back(NewSeq);
		SubsetMSA->m_Owners.push_back(true);
		}

	uint NewSeqCount = SIZE(SeqIndexes);
	vector<char> NewCol(NewSeqCount);
	for (uint OldColIndex = 0; OldColIndex < OldColCount; ++OldColIndex)
		{
		bool AllGaps = true;
		uint NewSeqIndex = 0;
		for (set<int>::const_iterator p = SeqIndexes.begin();
		  p != SeqIndexes.end(); ++p)
			{
			int OldSeqIndex = *p;
			const Sequence *OldSeq = GetSequence(OldSeqIndex);
			char c = OldSeq->GetChar(OldColIndex);
			NewCol[NewSeqIndex++] = c;
			if (c != '-')
				AllGaps = false;
			}
		if (AllGaps)
			continue;

		for (uint NewSeqIndex = 0; NewSeqIndex < NewSeqCount; ++NewSeqIndex)
			{
			char c = NewCol[NewSeqIndex];
			Sequence *NewSeq = (Sequence *) SubsetMSA->m_Seqs[NewSeqIndex];
			NewSeq->m_CharVec.push_back(c);
			}
		}

	return SubsetMSA;
	}
