#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "pwpath.h"
#include "textfile.h"
#include "timing.h"

SCORE AlignTwoMSAs(const MSA &msa1, const MSA &msa2, MSA &msaOut, PWPath &Path,
  bool bLockLeft, bool bLockRight)
	{
	const unsigned uLengthA = msa1.GetColCount();
	const unsigned uLengthB = msa2.GetColCount();

	ProfPos *PA = ProfileFromMSA(msa1);
	ProfPos *PB = ProfileFromMSA(msa2);

	if (bLockLeft)
		{
		PA[0].m_scoreGapOpen = MINUS_INFINITY;
		PB[0].m_scoreGapOpen = MINUS_INFINITY;
		}

	if (bLockRight)
		{
		PA[uLengthA-1].m_scoreGapClose = MINUS_INFINITY;
		PB[uLengthB-1].m_scoreGapClose = MINUS_INFINITY;
		}

	SCORE Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);

	AlignTwoMSAsGivenPath(Path, msa1, msa2, msaOut);

	delete[] PA;
	delete[] PB;

	return Score;
	}
