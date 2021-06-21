#ifndef ScoreHistory_h
#define ScoreHistory_h

class ScoreHistory
	{
public:
	ScoreHistory(unsigned uIters, unsigned uInternalNodeCount);
	~ScoreHistory();
	bool SetScore(unsigned uIter, unsigned uInternalNodeIndex, bool bRight, SCORE Score);
	void LogMe() const;
	SCORE GetScore(unsigned uIter, unsigned uInternalNodeIndex, bool bReversed,
	  bool bRight) const;

private:
	SCORE **m_Score;
	bool **m_bScoreSet;
	unsigned m_uIters;
	unsigned m_uNodeCount;
	};

#endif	// ScoreHistory_h
