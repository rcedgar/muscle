#pragma once

class Derep
	{
public:
	MultiSequence *m_InputSeqs = 0;
	vector<uint> m_SeqIndexToRepSeqIndex;
	vector<uint> m_RepSeqIndexes;
	vector<vector<uint> > m_RepSeqIndexToSeqIndexes;
	
	uint m_SlotCount = 0;
	vector<vector<uint> > m_HashToSeqIndexes;

	bool m_Disable = false;

public:
	uint CalcHash(const Sequence *Seq) const;
	void Clear();
	void Run(MultiSequence &InputSeqs, bool ShowProgress = true);
	uint Search(uint SeqIndex) const;
	void GetUniqueSeqs(MultiSequence &UniqueSeqs);
	void AddToHash(uint SeqIndex);
	bool SeqsEq(uint SeqIndex1, uint SeqIndex2) const;
	void Validate() const;
	void GetDupeGSIs(vector<uint> &GSIs,
	  vector<uint> &GlobalRepSeqIndexes) const;
	void GetRepLabelToDupeLabels(
	  unordered_map<string, vector<string> > &RepLabelToDupeLabels) const;
	};
