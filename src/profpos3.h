#pragma once

class MultiSequence;

class ProfPos3
	{
public:
	bool m_AllGaps;
	byte m_SortOrder[21];

// Frequencies
	float m_Freqs[20];

// Dimers are previous col + this col
	float m_LL;
	float m_LG;
	float m_GL;
	float m_GG;
	float m_fOcc;

// Scores
	float m_AAScores[20];
	float m_GapOpenScore;
	float m_GapCloseScore;

public:
	static uint m_NewCount;
	static uint m_DeleteCount;

public:
	ProfPos3()
		{
#pragma omp critical
		++m_NewCount;
		}

	~ProfPos3()
		{
#pragma omp critical
		++m_DeleteCount;
		}

	void SetFreqs(const MultiSequence &MSA, uint ColIndex,
	  const vector<float> &SeqWeights);
	void SetFreqs2(uint SeqCount,
	  uint LLCount, uint LGCount, uint GLCount, uint GGCount,
	  const vector<uint> &LetterCounts);
	void SetAAScores(const Mx2020 &SubstMxLetter);
	void SetOcc();
	void ToTsv(FILE *f) const;
	void LogMe() const;
	void SetStartDimers();
	//void NormalizeAAFreqsIfRequired();
	};
