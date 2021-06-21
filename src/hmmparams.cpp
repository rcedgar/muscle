#include "muscle.h"
#include "hmmparams.h"

const char *HMMTRANSToStr(HMMTRANS t)
	{
	switch (t)
		{
#define T(x) case HMMTRANS_##x: return #x;
#include "hmmtrans.h"
	default:
		asserta(false);
		}
	return "?";
	}

void HMMParams::GetProbs(HMMParams &Probs) const
	{
	if (m_Logs)
		ProbsToScores(*this, Probs);
	else
		Probs = *this;
	Probs.AssertProbsValid();
	}

void HMMParams::GetScores(HMMParams &Scores) const
	{
	if (m_Logs)
		Scores = *this;
	else
		{
		AssertProbsValid();
		ProbsToScores(*this, Scores);
		}
	}

void HMMParams::FromProbs(const HMMParams &Probs)
	{
	*this = Probs;
	AssertProbsValid();
	}

void HMMParams::FromScores(const HMMParams &Scores)
	{
	ScoresToProbs(Scores, *this);
	AssertProbsValid();
	}

void HMMParams::ScoresToProbs(const HMMParams &Scores, HMMParams &Probs)
	{
	for (uint i = 0; i < HMMTRANS_N; ++i)
		Probs.m_Trans[i] = EXP(Scores.m_Trans[i]);

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			Probs.m_Emits[i][j] = EXP(Scores.m_Emits[i][j]);

	Probs.m_Logs = false;
	Probs.AssertProbsValid();
	}

void HMMParams::ProbsToScores(const HMMParams &Probs, HMMParams &Scores)
	{
	Probs.AssertProbsValid();

	for (uint i = 0; i < HMMTRANS_N; ++i)
		Scores.m_Trans[i] = LOG(Probs.m_Trans[i]);

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			Scores.m_Emits[i][j] = LOG(Probs.m_Emits[i][j]);

	Scores.m_Logs = true;
	}

void HMMParams::AssertProbsValid() const
	{
	asserta(!m_Logs);
	asserta(strlen(HMM_ALPHA) == HMM_ALPHASIZE);
	asserta(SIZE(m_Trans) == HMMTRANS_N);
	asserta(SIZE(m_Emits) == HMM_ALPHASIZE);
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		asserta(SIZE(m_Emits[i]) == HMM_ALPHASIZE);

	float SumSTART = m_Trans[HMMTRANS_START_M] +
	  2*m_Trans[HMMTRANS_START_IS] +
	  2*m_Trans[HMMTRANS_START_IL];
	asserta(feq(SumSTART, 1.0));

	float SumM = m_Trans[HMMTRANS_M_M] +
	  2*m_Trans[HMMTRANS_M_IS] +
	  2*m_Trans[HMMTRANS_M_IL];
	asserta(feq(SumM, 1.0));

	float SumEmit = 0;
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			SumEmit += m_Emits[i][j];
	asserta(feq(SumEmit, 1.0));

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		for (uint j = 0; j < i; ++j)
			asserta(feq(m_Emits[i][j], m_Emits[j][i]));
	}

void HMMParams::ToFile(const string &FileName) const
	{
	if (FileName.empty())
		return;

	AssertProbsValid();
	FILE *f = CreateStdioFile(FileName);

#define T(x)	fprintf(f,	"T.%s	%.5g\n", #x, m_Trans[HMMTRANS_##x]);
#include "hmmtrans.h"

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char a = HMM_ALPHA[i];
		for (uint j = 0; j <= i; ++j)
			{
			char b = HMM_ALPHA[j];
			float P = m_Emits[i][j];
			fprintf(f, "E.%c%c	%.5g\n", a, b, P);
			}
		}

	CloseStdioFile(f);
	}

float HMMParams::GetNextProb(const string &Name)
	{
	string Line;
	vector<string> Fields;
	bool Ok = GetNextLine(Line);
	if (!Ok)
		Die("GetNextProb(%s) end-of-data", Name.c_str());
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2)
		Die("GetNextProb(%s) expected 2 fields got '%s'", Name.c_str(), Line.c_str());
	if (Fields[0] != Name)
		Die("ReadGetNextProbTrans(%s) got '%s'", Name.c_str(), Line.c_str());
	float P = (float) StrToFloat(Fields[1]);
	return P;
	}

bool HMMParams::GetNextLine(string &Line)
	{
	Line.clear();
	if (m_LineNr >= SIZE(m_Lines))
		return false;
	Line = m_Lines[m_LineNr++];
	return true;
	}

void HMMParams::FromFile(const string &FileName)
	{
	vector<string> Lines;
	ReadStringsFromFile(FileName, Lines);
	FromStrings(Lines);
	}

void HMMParams::FromStrings(const vector<string> &Lines)
	{
	m_Lines = Lines;
	m_LineNr = 0;

#define T(x)	m_Trans[HMMTRANS_##x] = GetNextProb("T." #x);
#include "hmmtrans.h"

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char a = HMM_ALPHA[i];
		for (uint j = 0; j <= i; ++j)
			{
			char b = HMM_ALPHA[j];

			string Name;
			Ps(Name, "E.%c%c", a, b);
			float P = GetNextProb(Name);
			m_Emits[i][j] = P;
			m_Emits[j][i] = P;
			}
		}

	AssertProbsValid();
	}

void HMMParams::FromDefaults()
	{
	vector<string> Lines;
	GetDefaultHMMParams(Lines);
	FromStrings(Lines);
	}

void HMMParams::ToPairHMM() const
	{
	HMMParams Scores;
	GetScores(Scores);

	HMMParams Probs;
	GetProbs(Probs);

	const vector<float> &Trans = Scores.m_Trans;
	const vector<vector<float> > &Emits = Scores.m_Emits;

	float SumInserts = 0;
	vector<float> InsertScores;
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		float MarginalProb = 0;
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			{
			float P = Probs.m_Emits[i][j];
			MarginalProb += P;
			}

		float Score = LOG(MarginalProb);
		InsertScores.push_back(Score);

		SumInserts += MarginalProb;
		}
	asserta(feq(SumInserts, 1.0));

	PairHMM::m_StartScore[HMMSTATE_M] = Trans[HMMTRANS_START_M];

	PairHMM::m_StartScore[HMMSTATE_ISX] = Trans[HMMTRANS_START_IS];
	PairHMM::m_StartScore[HMMSTATE_ISY] = Trans[HMMTRANS_START_IS];

	PairHMM::m_StartScore[HMMSTATE_ILX] = Trans[HMMTRANS_START_IL];
	PairHMM::m_StartScore[HMMSTATE_ILY] = Trans[HMMTRANS_START_IL];

	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		for (uint j = 0; j < HMMSTATE_COUNT; ++j)
			PairHMM::m_TransScore[i][j] = LOG_ZERO;

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_M] = Trans[HMMTRANS_M_M];

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_ISX] = Trans[HMMTRANS_M_IS];
	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_ISY] = Trans[HMMTRANS_M_IS];

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_ILX] = Trans[HMMTRANS_M_IL];
	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_ILY] = Trans[HMMTRANS_M_IL];

	PairHMM::m_TransScore[HMMSTATE_ISX][HMMSTATE_ISX] = Trans[HMMTRANS_IS_IS];
	PairHMM::m_TransScore[HMMSTATE_ISY][HMMSTATE_ISY] = Trans[HMMTRANS_IS_IS];

	PairHMM::m_TransScore[HMMSTATE_ILX][HMMSTATE_ILX] = Trans[HMMTRANS_IL_IL];
	PairHMM::m_TransScore[HMMSTATE_ILY][HMMSTATE_ILY] = Trans[HMMTRANS_IL_IL];

	PairHMM::m_TransScore[HMMSTATE_ISX][HMMSTATE_M] = Trans[HMMTRANS_IS_M];
	PairHMM::m_TransScore[HMMSTATE_ISY][HMMSTATE_M] = Trans[HMMTRANS_IS_M];

	PairHMM::m_TransScore[HMMSTATE_ILX][HMMSTATE_M] = Trans[HMMTRANS_IL_M];
	PairHMM::m_TransScore[HMMSTATE_ILY][HMMSTATE_M] = Trans[HMMTRANS_IL_M];

	for (uint i = 0; i < 256; ++i)
		PairHMM::m_InsScore[i] = LOG(WILDCARD_PROB);

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char a = HMM_ALPHA[i];
		byte ia = (byte) tolower(a);
		byte iA = (byte) toupper(a);
		float P = InsertScores[i];
		PairHMM::m_InsScore[ia] = P;
		PairHMM::m_InsScore[iA] = P;
		}

	for (uint i = 0; i < 256; ++i)
		for (uint j = 0; j < 256; ++j)
			PairHMM::m_MatchScore[i][j] = LOG(WILDCARD_PROB*WILDCARD_PROB);

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char a = HMM_ALPHA[i];
		byte ia = (byte) tolower(a);
		byte iA = (byte) toupper(a);
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			{
			float P = Emits[i][j];

			char b = HMM_ALPHA[j];
			byte ib = (byte) tolower(b);
			byte iB = (byte) toupper(b);

			PairHMM::m_MatchScore[ia][ib] = P;
			PairHMM::m_MatchScore[ia][iB] = P;
			PairHMM::m_MatchScore[iA][ib] = P;
			PairHMM::m_MatchScore[iA][iB] = P;
			}
		}
	}

void HMMParams::DeltaEmitProbs(const vector<vector<float> > &Deltas)
	{
	asserta(!m_Logs);
	asserta(SIZE(Deltas) == HMM_ALPHASIZE);

	float Sum = 0;
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		for (uint j = 0; j <= i; ++j)
			{
			float d = Deltas[i][j];
			asserta(d > 0);
			float P = m_Emits[i][j]*d;
			m_Emits[i][j] = P;
			m_Emits[j][i] = P;
			Sum += P;
			if (i != j)
				Sum += P;
			}
		}

	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			m_Emits[i][j] /= Sum;

	AssertProbsValid();
	}

void HMMParams::DeltaTransProbs(float dStartIS, float dStartIL,
  float dShortOpen, float dShortExtend,
  float dLongOpen, float dLongExtend)
	{
	DeltaStartProbs(dStartIS, dStartIL);
	DeltaShortGapProbs(dShortOpen, dShortExtend);
	DeltaLongGapProbs(dLongOpen, dLongExtend);
	}

void HMMParams::DeltaStartProbs(float dIS, float dIL)
	{
	asserta(!m_Logs);
	asserta(dIS >= 0);
	asserta(dIL >= 0);

	m_Trans[HMMTRANS_START_IS] *= dIS;
	m_Trans[HMMTRANS_START_IL] *= dIL;

	float Sum = 
	  m_Trans[HMMTRANS_START_M] +
	  2*m_Trans[HMMTRANS_START_IS] +
	  2*m_Trans[HMMTRANS_START_IL];

	m_Trans[HMMTRANS_START_M] /= Sum;
	m_Trans[HMMTRANS_START_IS] /= Sum;
	m_Trans[HMMTRANS_START_IL] /= Sum;

	AssertProbsValid();
	}

void HMMParams::DeltaShortGapProbs(float dOpen, float dExtend)
	{
	asserta(!m_Logs);
	asserta(dOpen > 0);
	asserta(dExtend >= 0);

	m_Trans[HMMTRANS_M_IS] *= dOpen;
	m_Trans[HMMTRANS_IS_IS] *= dExtend;

	float SumM = 
	  m_Trans[HMMTRANS_M_M] +
	  2*m_Trans[HMMTRANS_M_IS] +
	  2*m_Trans[HMMTRANS_M_IL];

	m_Trans[HMMTRANS_M_M] /= SumM;
	m_Trans[HMMTRANS_M_IS] /= SumM;
	m_Trans[HMMTRANS_M_IL] /= SumM;

	float SumIS =
	  m_Trans[HMMTRANS_IS_IS] +
	  m_Trans[HMMTRANS_IS_M];

	m_Trans[HMMTRANS_IS_IS] /= SumIS;
	m_Trans[HMMTRANS_IS_M] /= SumIS;

	AssertProbsValid();
	}

void HMMParams::DeltaLongGapProbs(float dOpen, float dExtend)
	{
	asserta(!m_Logs);
	asserta(dOpen > 0);
	asserta(dExtend >= 0);

	m_Trans[HMMTRANS_M_IL] *= dOpen;
	m_Trans[HMMTRANS_IL_IL] *= dExtend;

	float SumM = 
	  m_Trans[HMMTRANS_M_M] +
	  2*m_Trans[HMMTRANS_M_IS] +
	  2*m_Trans[HMMTRANS_M_IL];

	m_Trans[HMMTRANS_M_M] /= SumM;
	m_Trans[HMMTRANS_M_IS] /= SumM;
	m_Trans[HMMTRANS_M_IL] /= SumM;

	float SumIL =
	  m_Trans[HMMTRANS_IL_IL] +
	  m_Trans[HMMTRANS_IL_M];

	m_Trans[HMMTRANS_IL_IL] /= SumIL;
	m_Trans[HMMTRANS_IL_M] /= SumIL;

	AssertProbsValid();
	}
