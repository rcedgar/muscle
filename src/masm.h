#pragma once

#include "masmcol.h"

class MASM
	{
public:
	const MultiSequence *m_Aln = 0;
	string m_Label;
	uint m_ColCount = 0;
	uint m_SeqCount = 0;
	uint m_FeatureCount = 0;
	vector<string> m_FeatureNames;
	vector<uint> m_AlphaSizes;
	vector<MASMCol *> m_Cols;
	vector<string> m_UngappedSeqs;
	vector<vector<vector<byte> > > m_FeatureAlnVec;
	uint m_AAFeatureIdx = UINT_MAX;

public:
	void Clear()
		{
		m_Aln = 0;
		m_Label.clear();
		m_ColCount = 0;
		m_SeqCount = 0;
		m_FeatureNames.clear();
		m_AlphaSizes.clear();
		for (uint i = 0; i < SIZE(m_Cols); ++i)
			delete m_Cols[i];
		m_Cols.clear();
		m_UngappedSeqs.clear();
		m_FeatureAlnVec.clear();
		}

	uint GetColCount() const { return m_ColCount; }
	void FromMSA(const MultiSequence &Aln, const string &Label,
	  float GapOpen, float GapExt);
	const MASMCol &GetCol(uint i) const
		{
		assert(i < SIZE(m_Cols));
		return *m_Cols[i];
		}
	void SetUngappedSeqs();
	void SetFeatureAlnVec();
	void SetFeatureAln(uint FeatureIdx);
	void MakeSMx(const vector<vector<byte> > &ProfB, Mx<float> &SMx) const;
	void GetCounts(uint ColIndex, uint &LetterCount,
	  uint &GapOpenCount, uint &GapExtCount, uint &GapCloseCount);
	void GetFreqsVec(uint ColIndex, vector<vector<float> > &FreqsVec);
	void GetFreqs(uint ColIndex, uint FeatureIdx, vector<float> &Freqs);
	void ToFile(const string &FileName) const;
	void ToFile(FILE *f) const;
	void FromFile(const string &FileName);
	void MakeSMx_Sequence(const Sequence &Q, Mx<float> &SMx) const;
	void GetConsensusSeq(string &Seq) const;
	};
