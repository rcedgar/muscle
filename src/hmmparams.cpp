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

void HMMParams::FromParams(const HMMParams &Params, bool AsProbs)
	{
	if (AsProbs)
		{
		HMMParams Probs;
		Params.GetProbs(Probs);
		*this = Probs;
		}
	else
		{
		HMMParams Scores;
		Params.GetProbs(Scores);
		*this = Scores;
		}
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

void HMMParams::ToSingleAffineProbs(HMMParams &Params)
	{
	GetProbs(Params);

	vector<float> &T = Params.m_Trans;

	float SI = (T[HMMTRANS_START_IS] + T[HMMTRANS_START_IL])/2;

	float MI = (T[HMMTRANS_M_IS] + T[HMMTRANS_M_IL])/2;
	float IM = (T[HMMTRANS_IS_M] + T[HMMTRANS_IL_M])/2;
	float II = (T[HMMTRANS_IS_IS] + T[HMMTRANS_IL_IL])/2;

	T[HMMTRANS_START_IS] = SI;
	T[HMMTRANS_START_IL] = SI;

	T[HMMTRANS_M_IS] = MI;
	T[HMMTRANS_M_IL] = MI;

	T[HMMTRANS_IS_M] = IM;
	T[HMMTRANS_IL_M] = IM;

	T[HMMTRANS_IS_IS] = II;
	T[HMMTRANS_IL_IL] = II;

	Params.AssertProbsValid();
	}

void HMMParams::ScoresToProbs(const HMMParams &Scores, HMMParams &Probs)
	{
	Probs.m_Alpha = Scores.m_Alpha;
	const unsigned AlphaSize = Scores.GetAlphaSize();

	for (uint i = 0; i < HMMTRANS_N; ++i)
		Probs.m_Trans[i] = exp(Scores.m_Trans[i]);

	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j < AlphaSize; ++j)
			Probs.m_Emits[i][j] = exp(Scores.m_Emits[i][j]);

	Probs.m_Logs = false;
	Probs.AssertProbsValid();
	}

void HMMParams::ProbsToScores(const HMMParams &Probs, HMMParams &Scores)
	{
	Probs.AssertProbsValid();
	Scores.m_Alpha = Probs.m_Alpha;
	const unsigned AlphaSize = Probs.GetAlphaSize();

	Scores.m_Trans.clear();
	Scores.m_Trans.resize(HMMTRANS_N);

	Scores.m_Emits.clear();
	Scores.m_Emits.clear();
	Scores.m_Emits.resize(AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		Scores.m_Emits[i].resize(AlphaSize, FLT_MAX);

	for (uint i = 0; i < HMMTRANS_N; ++i)
		Scores.m_Trans[i] = log(Probs.m_Trans[i]);

	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j < AlphaSize; ++j)
			Scores.m_Emits[i][j] = log(Probs.m_Emits[i][j]);

	Scores.m_Logs = true;
	}

void HMMParams::AssertProbsValid() const
	{
	asserta(!m_Logs);
	asserta(SIZE(m_Trans) == HMMTRANS_N);
	const uint AlphaSize = GetAlphaSize();
	asserta(SIZE(m_Emits) == AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		asserta(SIZE(m_Emits[i]) == AlphaSize);

	float SumSTART = m_Trans[HMMTRANS_START_M] +
	  2*m_Trans[HMMTRANS_START_IS] +
	  2*m_Trans[HMMTRANS_START_IL];
	asserta(feq(SumSTART, 1.0));

	float SumIS = m_Trans[HMMTRANS_IS_M] + m_Trans[HMMTRANS_IS_IS];
	asserta(feq(SumIS, 1.0));

	float SumIL = m_Trans[HMMTRANS_IL_M] + m_Trans[HMMTRANS_IL_IL];
	asserta(feq(SumIL, 1.0));

	float SumM = m_Trans[HMMTRANS_M_M] +
	  2*m_Trans[HMMTRANS_M_IS] +
	  2*m_Trans[HMMTRANS_M_IL];
	asserta(feq(SumM, 1.0));

	float SumEmit = 0;
	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j < AlphaSize; ++j)
			SumEmit += m_Emits[i][j];
	asserta(feq(SumEmit, 1.0));

	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j < i; ++j)
			asserta(feq(m_Emits[i][j], m_Emits[j][i]));
	}

void HMMParams::ToFile(const string &FileName) const
	{
	if (FileName.empty())
		return;

	const uint AlphaSize = GetAlphaSize();
	AssertProbsValid();
	FILE *f = CreateStdioFile(FileName);

	if (m_Alpha == AMINO_ALPHA)
		fprintf(f, "HMM	aa\n");
	else if (m_Alpha == NT_ALPHA)
		fprintf(f, "HMM	nt\n");
	else
		Die("HMMParams::ToFile alpha='%'s", m_Alpha.c_str());

#define T(x)	fprintf(f,	"T.%s	%.5g\n", #x, m_Trans[HMMTRANS_##x]);
#include "hmmtrans.h"

	for (uint i = 0; i < AlphaSize; ++i)
		{
		char a = m_Alpha[i];
		for (uint j = 0; j <= i; ++j)
			{
			char b = m_Alpha[j];
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

	string AlphaLine;
	bool Ok = GetNextLine(AlphaLine);
	asserta(Ok);

	vector<string> Fields;
	Split(AlphaLine, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "HMM")
		Die("Invalid HMM file");
	m_Alpha = Fields[1];
	if (Fields[1] == "aa")
		m_Alpha = AMINO_ALPHA;
	else if (Fields[1] == "nt")
		m_Alpha = NT_ALPHA;
	else
		Die("Invalid HMM alphabet '%s'", m_Alpha.c_str());

	m_Trans.clear();
	m_Trans.resize(HMMTRANS_N, FLT_MAX);

#define T(x)	m_Trans[HMMTRANS_##x] = GetNextProb("T." #x);
#include "hmmtrans.h"

	const uint AlphaSize = SIZE(m_Alpha);
	m_Emits.clear();
	m_Emits.resize(AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		m_Emits[i].resize(AlphaSize, FLT_MAX);

	for (uint i = 0; i < AlphaSize; ++i)
		{
		char a = m_Alpha[i];
		for (uint j = 0; j <= i; ++j)
			{
			char b = m_Alpha[j];

			string Name;
			Ps(Name, "E.%c%c", a, b);
			float P = GetNextProb(Name);
			m_Emits[i][j] = P;
			m_Emits[j][i] = P;
			}
		}

	Normalize();
	AssertProbsValid();
	}

void HMMParams::FromDefaults(bool Nucleo)
	{
	vector<string> Lines;
	if (Nucleo)
		GetDefaultHMMParams_Nucleo(Lines);
	else
		GetDefaultHMMParams_Amino(Lines);
	FromStrings(Lines);
	}

void HMMParams::CmdLineUpdate()
	{
	if (!optset_m_is && !optset_m_il && !optset_is_is && !optset_il_il
	  && !optset_s_is && !optset_s_il)
		return;
	asserta(!m_Logs);
	if (optset_s_is)	m_Trans[HMMTRANS_START_IS] = (float) opt(s_is);
	if (optset_s_il)	m_Trans[HMMTRANS_START_IL] = (float) opt(s_il);
	if (optset_m_is)	m_Trans[HMMTRANS_M_IS] = (float) opt(m_is);
	if (optset_m_il)	m_Trans[HMMTRANS_M_IL] = (float) opt(m_il);
	if (optset_is_is)	m_Trans[HMMTRANS_IS_IS] = (float) opt(is_is);
	if (optset_il_il)	m_Trans[HMMTRANS_IL_IL] = (float) opt(il_il);
	Normalize();
	}

void HMMParams::ToPairHMM() const
	{
	HMMParams Scores;
	GetScores(Scores);

	HMMParams Probs;
	GetProbs(Probs);

	const uint AlphaSize = GetAlphaSize();

	const vector<float> &Trans = Scores.m_Trans;
	const vector<vector<float> > &Emits = Scores.m_Emits;

	float SumInserts = 0;
	vector<float> InsertScores;
	for (uint i = 0; i < AlphaSize; ++i)
		{
		float MarginalProb = 0;
		for (uint j = 0; j < AlphaSize; ++j)
			{
			float P = Probs.m_Emits[i][j];
			MarginalProb += P;
			}

		float Score = log(MarginalProb);
		InsertScores.push_back(Score);

		SumInserts += MarginalProb;
		}
	asserta(feq(SumInserts, 1.0));

	PairHMM::m_StartScore[HMMSTATE_M] = Trans[HMMTRANS_START_M];

	PairHMM::m_StartScore[HMMSTATE_IX] = Trans[HMMTRANS_START_IS];
	PairHMM::m_StartScore[HMMSTATE_IY] = Trans[HMMTRANS_START_IS];

	PairHMM::m_StartScore[HMMSTATE_JX] = Trans[HMMTRANS_START_IL];
	PairHMM::m_StartScore[HMMSTATE_JY] = Trans[HMMTRANS_START_IL];

	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		for (uint j = 0; j < HMMSTATE_COUNT; ++j)
			PairHMM::m_TransScore[i][j] = LOG_ZERO;

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_M] = Trans[HMMTRANS_M_M];

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_IX] = Trans[HMMTRANS_M_IS];
	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_IY] = Trans[HMMTRANS_M_IS];

	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_JX] = Trans[HMMTRANS_M_IL];
	PairHMM::m_TransScore[HMMSTATE_M][HMMSTATE_JY] = Trans[HMMTRANS_M_IL];

	PairHMM::m_TransScore[HMMSTATE_IX][HMMSTATE_IX] = Trans[HMMTRANS_IS_IS];
	PairHMM::m_TransScore[HMMSTATE_IY][HMMSTATE_IY] = Trans[HMMTRANS_IS_IS];

	PairHMM::m_TransScore[HMMSTATE_JX][HMMSTATE_JX] = Trans[HMMTRANS_IL_IL];
	PairHMM::m_TransScore[HMMSTATE_JY][HMMSTATE_JY] = Trans[HMMTRANS_IL_IL];

	PairHMM::m_TransScore[HMMSTATE_IX][HMMSTATE_M] = Trans[HMMTRANS_IS_M];
	PairHMM::m_TransScore[HMMSTATE_IY][HMMSTATE_M] = Trans[HMMTRANS_IS_M];

	PairHMM::m_TransScore[HMMSTATE_JX][HMMSTATE_M] = Trans[HMMTRANS_IL_M];
	PairHMM::m_TransScore[HMMSTATE_JY][HMMSTATE_M] = Trans[HMMTRANS_IL_M];

	float WildcardInsertProb = 1.0f/AlphaSize;
	for (uint i = 0; i < 256; ++i)
		PairHMM::m_InsScore[i] = log(WildcardInsertProb);

	for (uint i = 0; i < AlphaSize; ++i)
		{
		char a = m_Alpha[i];
		byte ia = (byte) tolower(a);
		byte iA = (byte) toupper(a);
		float P = InsertScores[i];
		PairHMM::m_InsScore[ia] = P;
		PairHMM::m_InsScore[iA] = P;
		}

	for (uint i = 0; i < 256; ++i)
		for (uint j = 0; j < 256; ++j)
			PairHMM::m_MatchScore[i][j] = log(WildcardInsertProb*WildcardInsertProb);

	for (uint i = 0; i < AlphaSize; ++i)
		{
		char a = m_Alpha[i];
		byte ia = (byte) tolower(a);
		byte iA = (byte) toupper(a);
		for (uint j = 0; j < AlphaSize; ++j)
			{
			float P = Emits[i][j];

			char b = m_Alpha[j];
			byte ib = (byte) tolower(b);
			byte iB = (byte) toupper(b);

			PairHMM::m_MatchScore[ia][ib] = P;
			PairHMM::m_MatchScore[ia][iB] = P;
			PairHMM::m_MatchScore[iA][ib] = P;
			PairHMM::m_MatchScore[iA][iB] = P;
			}
		}

	if (optset_anchor_letter)
		{
		const string sAL = opt(anchor_letter);
		asserta(SIZE(sAL) == 1);
		char c = sAL[0];
		PairHMM::m_MatchScore[c][c] = -0.1f;
		}

	if (AlphaSize == 4)
		PairHMM::FixUT();
	}

void PairHMM::FixUT()
	{
	PairHMM::m_InsScore['U'] = PairHMM::m_InsScore['T'];
	PairHMM::m_InsScore['u'] = PairHMM::m_InsScore['t'];

	for (uint i = 0; i < 256; ++i)
		{
		float P = PairHMM::m_MatchScore['T'][i];

		PairHMM::m_MatchScore['U'][i] = P;
		PairHMM::m_MatchScore['u'][i] = P;
		PairHMM::m_MatchScore[i]['U'] = P;
		PairHMM::m_MatchScore[i]['u'] = P;
		}
	}

void HMMParams::NormalizeEmit()
	{
	asserta(!m_Logs);

	const uint AlphaSize = GetAlphaSize();
	float Sum = 0;
	for (uint i = 0; i < AlphaSize; ++i)
		{
		for (uint j = 0; j <= i; ++j)
			{
			float P = m_Emits[i][j];
			m_Emits[i][j] = P;
			m_Emits[j][i] = P;
			Sum += P;
			if (i != j)
				Sum += P;
			}
		}

	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j < AlphaSize; ++j)
			m_Emits[i][j] /= Sum;
	}

void HMMParams::NormalizeStart()
	{
	float Sum = 
	  m_Trans[HMMTRANS_START_M] +
	  2*m_Trans[HMMTRANS_START_IS] +
	  2*m_Trans[HMMTRANS_START_IL];

	m_Trans[HMMTRANS_START_M] /= Sum;
	m_Trans[HMMTRANS_START_IS] /= Sum;
	m_Trans[HMMTRANS_START_IL] /= Sum;
	}

void HMMParams::NormalizeShortGap()
	{
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
	}

void HMMParams::NormalizeLongGap()
	{
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
	}

void HMMParams::NormalizeMatch()
	{
	float SumM = m_Trans[HMMTRANS_M_M] +
	  2*m_Trans[HMMTRANS_M_IS] +
	  2*m_Trans[HMMTRANS_M_IL];
	
	m_Trans[HMMTRANS_M_M] /= SumM;
	m_Trans[HMMTRANS_M_IS] /= SumM;
	m_Trans[HMMTRANS_M_IL] /= SumM;
	}

void HMMParams::Normalize()
	{
	NormalizeStart();
	NormalizeShortGap();
	NormalizeLongGap();
	NormalizeEmit();

	AssertProbsValid();
	}
