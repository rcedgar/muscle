#include "myutils.h"
#include "pairhmm.h"
#include "hmmparams.h"

float PairHMM::m_StartScore[HMMSTATE_COUNT];
float PairHMM::m_TransScore[HMMSTATE_COUNT][HMMSTATE_COUNT];
float PairHMM::m_MatchScore[256][256];
float PairHMM::m_InsScore[256];

/***
                    M        IX        IY        JX        JY
                  [0]       [1]       [2]       [3]       [4]
  M   [0]       0.960     0.012     0.012   0.00801   0.00801
 IX   [1]       0.603     0.397         0         0         0
 IY   [2]       0.603         0     0.397         0         0
 JX   [3]       0.101         0         0     0.899         0
 JY   [4]       0.101         0         0         0     0.899
***/
static vector<vector<float> > transMat;

static void ConstructTransMat(const vector<float>& gapOpen,
  const vector<float>& gapExtend)
	{
	transMat.clear();
	transMat.resize(HMMSTATE_COUNT);
	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		transMat[i].resize(HMMSTATE_COUNT);

	transMat[0][0] = 1;
	for (uint i = 0; i < InsertStateCount; i++)
		{
		transMat[0][2 * i + 1] = gapOpen[2 * i];
		transMat[0][2 * i + 2] = gapOpen[2 * i + 1];

		transMat[0][0] -= (gapOpen[2 * i] + gapOpen[2 * i + 1]);

		transMat[2 * i + 1][2 * i + 1] = gapExtend[2 * i];
		transMat[2 * i + 2][2 * i + 2] = gapExtend[2 * i + 1];

		transMat[2 * i + 1][2 * i + 2] = 0;
		transMat[2 * i + 2][2 * i + 1] = 0;

		transMat[2 * i + 1][0] = 1 - gapExtend[2 * i];
		transMat[2 * i + 2][0] = 1 - gapExtend[2 * i + 1];
		}
	asserta(transMat[0][0] > 0);
	}

void PairHMM::Create2(const vector<float>& initDistribMat,
  const vector<vector<float> > &transMat, const vector<float>& emitSingle,
  const vector<vector<float> > &emitPairs)
	{
	for (int i = 0; i < HMMSTATE_COUNT; i++)
		m_StartScore[i] = log(initDistribMat[i]);

	for (int i = 0; i < HMMSTATE_COUNT; i++)
		for (int j = 0; j < HMMSTATE_COUNT; j++)
			m_TransScore[i][j] = log(transMat[i][j]);

	for (int i = 0; i < 256; i++)
		m_InsScore[i] = log(emitSingle[i]);

	for (int i = 0; i < 256; i++)
		for (int j = 0; j < 256; j++)
			m_MatchScore[i][j] = log(emitPairs[i][j]);
	}

void PairHMM::Create(const vector<float>& initDistribMat,
  const vector<float>& gapOpen, const vector<float>& gapExtend,
  const vector<vector<float>>& emitPairs, const vector<float>& emitSingle)
	{
	ConstructTransMat(gapOpen,  gapExtend);
	Create2(initDistribMat, transMat, emitSingle, emitPairs);
	}
