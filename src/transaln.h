#pragma once

class TransAln
	{
public:
// Input data
	const MultiSequence *m_MSA = 0;
	const MultiSequence *m_FreshSeqs = 0;
	const vector<uint> *m_FreshIndexToMSAIndex = 0;
	const vector<string> *m_PWPaths = 0;
	uint m_MSAColCount = 0;

// Derived data
	uint m_ExtendedMSAColCount = 0;
	vector<Sequence *> m_UngappedMSASeqs;
	vector<string> m_MSAPaths;
	vector<string> m_TPaths1;
	vector<string> m_TPaths2;
	vector<uint> m_MSAColToMaxInserts;
	string m_MPath;
	MultiSequence *m_ExtendedMSA = 0;

public:
	void Init(const MultiSequence &MSA,
	  const MultiSequence &FreshSeqs,
	  const vector<uint> &FreshIndexToMSAIndex,
	  const vector<string> &PWPaths);

	void LogMe() const;
	void LogTPath1Aln(uint FreshIndex) const;
	void LogTPath2Aln(uint FreshIndex, bool WithPath = true) const;
	void LogMPathAln(uint MSAIndex, bool WithPath = true) const;

	uint GetFreshCount() const;
	uint GetMSACount() const;
	const string &GetFreshLabel(uint MSAIndex) const;
	const string &GetMSALabel(uint MSAIndex) const;
	const Sequence &GetMSASeq(uint MSAIndex) const;
	const Sequence &GetUngappedMSASeq(uint MSAIndex) const;
	const Sequence &GetFreshSeq(uint FreshIndex) const;
	uint GetMSAIndex(uint FreshIndex) const;
	const string &GetPWPath(uint FreshIndex) const;
	const string &GetMSAPath(uint FreshIndex) const;
	const string &GetTPath1(uint FreshIndex) const;
	void MakeMPath(string &MPath) const;
	const string &GetTPath2(uint FreshIndex) const;
	uint GetUngappedMSASeqLength(uint MSAIndex) const;
	uint GetFreshSeqLength(uint MSAIndex) const;
	uint GetMSAColCount() const { return m_MSAColCount; }
	void MakeMSAPath(uint MSAIndex, string &MSAPath) const;
	void MakeTPath1(uint FreshIndex, string &Path1) const;
	void MakeTPath2(uint FreshIndex, string &Path2) const;
	void MakeMSAColToInserts(uint FreshIndex, vector<uint> &MSAColToInserts) const;
	void SetMaxInserts();
	void MakeExtendedMSA();
	Sequence *ExtendMSASeq(uint MSAIndex) const;
	Sequence *ExtendFreshSeq(uint MSAIndex) const;
	};
