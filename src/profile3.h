#pragma once

#include "profpos3.h"
#include "multisequence.h"

class Profile3
	{
public:
	vector<ProfPos3 *> m_PPs;

public:
	static uint m_NewCount;
	static uint m_DeleteCount;

public:
	Profile3()
		{
#pragma omp critical
		++m_NewCount;
		}

	~Profile3()
		{
#pragma omp critical
		++m_DeleteCount;
		Clear();
		}

public:
	void Clear();
	uint GetColCount() const { return SIZE(m_PPs); }
	void FromMSA(const MultiSequence &MSA,
	  const Mx2020 &SubstMx_Letter, float GapOpen,
	  const vector<float> &SeqWeights);
	void FromSeq(const Sequence &Seq,
	  const Mx2020 &SubstMx_Letter, float GapOpen);
	const ProfPos3 &GetPP(uint ColIndex) const;
	void SetScores(const Mx2020 &SubstMx_Letter, float GapOpen);
	void LogMe(const MultiSequence *MSA = 0) const;
	void Validate() const;
	void ToTsv(FILE *f) const;
	void ToTsv(const string &FileName) const;
	void SetAAScores(const Mx2020 &SubstMx_Letter);
	void SetGapScores(float GapOpen);
	void SetGapOpenScore(float GapOpen, uint ColIndex);
	void SetGapCloseScore(float GapOpen, uint ColIndex);
	uint LogDiffs(const Profile3 &Prof2) const;
	float GetSelfScore() const;
	};
