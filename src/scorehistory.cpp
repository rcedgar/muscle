#include "muscle.h"
#include "scorehistory.h"
#include <stdio.h>

#define TRACE	0

ScoreHistory::ScoreHistory(unsigned uIters, unsigned uNodeCount)
	{
	m_uNodeCount = uNodeCount;
	m_uIters = uIters;

	m_Score = new SCORE *[uIters];
	m_bScoreSet = new bool *[uIters];
	for (unsigned n = 0; n < uIters; ++n)
		{
		m_Score[n] = new SCORE[uNodeCount*2];
		m_bScoreSet[n] = new bool[uNodeCount*2];
		memset(m_bScoreSet[n], 0, uNodeCount*2*sizeof(bool));
		}
	}

ScoreHistory::~ScoreHistory()
	{
	for (unsigned n = 0; n < m_uIters; ++n)
		{
		delete[] m_Score[n];
		delete[] m_bScoreSet[n];
		}
	delete[] m_Score;
	delete[] m_bScoreSet;
	}

bool ScoreHistory::SetScore(unsigned uIter, unsigned uNodeIndex, bool bRight, SCORE Score)
	{
#if	TRACE
	Log("ScoreHistory::SetScore(Iter=%u Node=%u Right=%d Score=%g)\n",
	  uIter, uNodeIndex, bRight, Score);
#endif
	if (uIter >= m_uIters)
		Quit("ScoreHistory::SetScore-1");
	if (uNodeIndex >= m_uNodeCount)
		Quit("ScoreHistory::SetScore-2");

	const unsigned uIndex = uNodeIndex*2 + bRight;
	for (unsigned n = 1; n < uIter; ++n)
		{
		const unsigned uPrevIter = n - 1;
		if (!m_bScoreSet[uPrevIter][uIndex])
			{
			LogMe();
			Quit("ScoreHistory::SetScore-3");
			}
		if (m_Score[uPrevIter][uIndex] == Score)
			{
			ProgressStepsDone();
#if	TRACE
			Log("Oscillating\n");
#endif
			return true;
			}
		}
	m_Score[uIter][uIndex] = Score;
	m_bScoreSet[uIter][uIndex] = true;
	return false;
	}

void ScoreHistory::LogMe() const
	{
	Log("ScoreHistory\n");
	Log("Iter  Node  Right      Score\n");
	Log("----  ----  -----  ---------\n");
	for (unsigned uIter = 0; uIter < m_uIters; ++uIter)
		{
		bool bAnySet = false;
		for (unsigned n = 0; n < m_uNodeCount*2; ++n)
			if (m_bScoreSet[uIter][n])
				{
				bAnySet = true;
				break;
				}
		if (!bAnySet)
			return;
		for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
			{
			const unsigned uBase = 2*uNodeIndex;
			if (m_bScoreSet[uIter][uBase])
				Log("%4u  %4u         F  %9.3f\n", uIter, uNodeIndex, m_Score[uIter][uBase]);
			if (m_bScoreSet[uIter][uBase+1])
				Log("%4u  %4u         T  %9.3f\n", uIter, uNodeIndex, m_Score[uIter][uBase+1]);
			}
		}
	}

SCORE ScoreHistory::GetScore(unsigned uIter, unsigned uNodeIndex,
  bool bReverse, bool bRight) const
	{
	const unsigned uIndex = uNodeIndex*2 + bRight;
	if (!m_bScoreSet[uIter][uIndex])
		Quit("ScoreHistory::GetScore");
	return m_Score[uIter][uIndex];
	}
