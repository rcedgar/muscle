#include "muscle.h"
#include "msa.h"
#include "objscore.h"
#include "profile.h"
#include "timing.h"

#if	TIMING
TICKS g_ticksObjScore = 0;
#endif

SCORE ObjScore(const MSA &msa, const unsigned SeqIndexes1[],
  unsigned uSeqCount1, const unsigned SeqIndexes2[], unsigned uSeqCount2)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	const unsigned uSeqCount = msa.GetSeqCount();

	OBJSCORE OS = g_ObjScore;
	if (g_ObjScore == OBJSCORE_SPM)
		{
        if (uSeqCount <= 100)
			OS = OBJSCORE_XP;
		else
			OS = OBJSCORE_SPF;
		}

	MSA msa1;
	MSA msa2;

	switch (OS)
		{
	case OBJSCORE_DP:
	case OBJSCORE_XP:
		MSAFromSeqSubset(msa, SeqIndexes1, uSeqCount1, msa1);
		MSAFromSeqSubset(msa, SeqIndexes2, uSeqCount2, msa2);

		SetMSAWeightsMuscle(msa1);
		SetMSAWeightsMuscle(msa2);
		break;

	case OBJSCORE_SP:
	case OBJSCORE_SPF:
	case OBJSCORE_PS:
	// Yuck -- casting away const (design flaw)
		SetMSAWeightsMuscle((MSA &) msa);
		break;
		}

	SCORE Score = 0;
	switch (OS)
		{
	case OBJSCORE_SP:
		Score = ObjScoreSP(msa);
		break;

	case OBJSCORE_DP:
		Score = ObjScoreDP(msa1, msa2);
		break;

	case OBJSCORE_XP:
		Score = ObjScoreXP(msa1, msa2);
		break;

	case OBJSCORE_PS:
		Score = ObjScorePS(msa);
		break;

	case OBJSCORE_SPF:
		Score = ObjScoreSPDimer(msa);
		break;
	
	default:
		Quit("Invalid g_ObjScore=%d", g_ObjScore);
		}
#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksObjScore += (t2 - t1);
#endif
	return Score;
	}

SCORE ObjScoreIds(const MSA &msa, const unsigned Ids1[],
  unsigned uCount1, const unsigned Ids2[], unsigned uCount2)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	unsigned *SeqIndexes1 = new unsigned[uCount1];
	unsigned *SeqIndexes2 = new unsigned[uCount2];

	for (unsigned n = 0; n < uCount1; ++n)
		SeqIndexes1[n] = msa.GetSeqIndex(Ids1[n]);

	for (unsigned n = 0; n < uCount2; ++n)
		SeqIndexes2[n] = msa.GetSeqIndex(Ids2[n]);

#if DOUBLE_AFFINE
	extern SCORE ObjScoreDA(const MSA &msa, SCORE *ptrLetters, SCORE *ptrGaps);
	SCORE Letters, Gaps;
	SCORE dObjScore = ObjScoreDA(msa, &Letters, &Gaps);

	delete[] SeqIndexes1;
	delete[] SeqIndexes2;
#else
	SCORE dObjScore = ObjScore(msa, SeqIndexes1, uCount1, SeqIndexes2, uCount2);
#endif
#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksObjScore += (t2 - t1);
#endif
	return dObjScore;
	}
