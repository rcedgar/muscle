#pragma once

class MASM;

class MASMCol
	{
public:
	const MASM *m_MASM = 0;

	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	float m_GapClose = FLT_MAX;

	float m_LetterFreq = FLT_MAX;
	float m_GapOpenFreq = FLT_MAX;
	float m_GapExtFreq = FLT_MAX;
	float m_GapCloseFreq = FLT_MAX;

// Freq = m_FreqsVec[FeatureIdx][Letter]
	vector<vector<float> > m_FreqsVec;

// Scores derived from Freqs which includes gaps,
//   so scores already account for occupancy
// Score = m_ScoresVec[FeatureIdx][Letter]
	vector<vector<float> > m_ScoresVec;

//// Letter = m_SortOrderVec[FeatureIdx][k]
//	vector<vector<byte> > m_SortOrderVec;

public:
	void SetScoreVec();
	//void SetSortOrderVec();
	void ToFile(FILE *f, uint ColIndex) const;
	const vector<float> &GetAAScores() const;
	};
