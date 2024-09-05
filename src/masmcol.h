#pragma once

class MASM;

class MASMCol
	{
public:
	const MASM *m_MASM = 0;
	uint m_ColIndex = 0;
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

public:
	void SetScoreVec();
	void ToFile(FILE *f) const;
	void FromFile(FILE *f);
	const vector<float> &GetAAScores() const;
	char GetConsensusAAChar() const;
	float GetMatchScore_MegaProfilePos(const vector<byte> &ProfPos) const;
	void LogMe() const;
	};
