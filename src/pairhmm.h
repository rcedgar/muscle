#ifndef pairhmm_h
#define pairhmm_h

#include "scoretype.h"
#include "multisequence.h"
#include "hmmparams.h"

static const uint DEFAULT_CONSISTENCY_ITERS = 2;
static const uint DEFAULT_REFINE_ITERS = 100;

enum HMMSTATE
	{
	HMMSTATE_M = 0,
	HMMSTATE_IX = 1,
	HMMSTATE_IY = 2,
	HMMSTATE_JX = 3,
	HMMSTATE_JY = 4,
	HMMSTATE_COUNT = 5
	};

static const uint InsertStateCount = 2;

class PairHMM
	{
public:
	static float m_StartScore[HMMSTATE_COUNT];
	static float m_TransScore[HMMSTATE_COUNT][HMMSTATE_COUNT];
	static float m_MatchScore[256][256];
	static float m_InsScore[256];

public:
	static void Create(const vector<float>& initDistribMat,
	  const vector<float>& gapOpen, const vector<float>& gapExtend,
	  const vector<vector<float>>& emitPairs,
	  const vector<float>& emitSingle);

	static void Create2(const vector<float>& initDistribMat,
	  const vector<vector<float> > &transMat, const vector<float>& emitSingle,
	  const vector<vector<float> > &emitPairs);

	static void WriteParamsReport(const string &FileName);
	static void WriteParamsReport(FILE *f);
	static void FixUT();
	};

#endif
