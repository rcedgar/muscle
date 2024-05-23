const float tSM = PairHMM::m_StartScore[HMMSTATE_M];
const float tSI = PairHMM::m_StartScore[HMMSTATE_IX];
const float tSJ = PairHMM::m_StartScore[HMMSTATE_JX];

const float tMM = PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_M];
const float tMI = PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_IX];
const float tMJ = PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_JX];

const float tII = PairHMM::m_TransScore[HMMSTATE_IX][HMMSTATE_IX];
const float tIM = PairHMM::m_TransScore[HMMSTATE_IX][HMMSTATE_M];

const float tJJ = PairHMM::m_TransScore[HMMSTATE_JX][HMMSTATE_JX];
const float tJM = PairHMM::m_TransScore[HMMSTATE_JX][HMMSTATE_M];

static const t_ByteMx &MatchScore = PairHMM::m_MatchScore;
static const t_ByteVec &InsScore = PairHMM::m_InsScore;
