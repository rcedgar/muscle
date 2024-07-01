#include "myutils.h"
#include "mega.h"
#include "alpha.h"
#include "pairhmm.h"

static const float VERY_SMALL_FREQ = 1e-6f;

const string &Mega::GetNextLine()
	{
	asserta(m_NextLineNr < SIZE(m_Lines));
	return m_Lines[m_NextLineNr++];
	}

void Mega::GetNextFields(vector<string> &Fields,
  uint ExpectedNrFields)
	{
	const string &Line = GetNextLine();
	Split(Line, Fields, '\t');
	if (ExpectedNrFields != UINT_MAX && SIZE(Fields) != ExpectedNrFields)
		Die("Expected %u fields got %u in '%s'",
		  ExpectedNrFields, SIZE(Fields), Line.c_str());
	}

void Mega::AssertSymmetrical(const vector<vector<float> > &Mx) const
	{
	const uint N = SIZE(Mx);
	for (uint i = 0; i < N; ++i)
		{
		const vector<float> &Row = Mx[i];
		asserta(SIZE(Row) == N);
		for (uint j = 0; j < i; ++j)
			asserta(feq(Mx[i][j], Mx[j][i]));
		}
	}

void Mega::CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
  vector<float> &MarginalFreqs) const
	{
	MarginalFreqs.clear();
	AssertSymmetrical(FreqsMx);
	const uint N = SIZE(FreqsMx);
	float SumMarginalFreqs = 0;
	for (uint i = 0; i < N; ++i)
		{
		const vector<float> &Row = FreqsMx[i];
		float MarginalFreq = 0;
		for (uint j = 0; j < N; ++j)
			MarginalFreq += Row[j];
		MarginalFreqs.push_back(MarginalFreq);
		SumMarginalFreqs += MarginalFreq;
		}
	asserta(feq(SumMarginalFreqs, 1));
	}

void Mega::CalcLogProbsMx(const vector<vector<float > > &FreqsMx,
  vector<vector<float > > &LogProbMx) const
	{
	AssertSymmetrical(FreqsMx);
	const uint N = SIZE(FreqsMx);
	LogProbMx.clear();
	LogProbMx.resize(N);
	for (uint i = 0; i < N; ++i)
		LogProbMx[i].resize(N);

	const float VERY_SMALL_FREQ = 1e-6f;
	for (uint i = 0; i < N; ++i)
		{
		const vector<float> &FreqsRow = FreqsMx[i];
		vector<float> &LogProbsRow = LogProbMx[i];
		for (uint j = 0; j < N; ++j)
			{
			float Freq = FreqsRow[j];
			asserta(Freq <= 1);
			if (Freq < VERY_SMALL_FREQ)
				Freq = VERY_SMALL_FREQ;
			LogProbsRow[j] = logf(Freq);
			}
		}
	AssertSymmetrical(LogProbMx);
	}

void Mega::FromFile(const string &FileName)
{
    FILE *f = OpenStdioFile(FileName);
    string Line;
    while (ReadLineStdioFile(f, Line))
        m_Lines.push_back(Line);
    CloseStdioFile(f);
    
    vector<string> flds;
    GetNextFields(flds, 3);
    asserta(flds[0] == "mega");
    m_FeatureCount = StrToUint(flds[1]);
    uint ProfileCount = StrToUint(flds[2]);
    m_LogProbsVec.resize(m_FeatureCount);
    m_LogProbMxVec.resize(m_FeatureCount);
    for (uint FeatureIdx = 0; FeatureIdx < m_FeatureCount; ++FeatureIdx)
    {
        GetNextFields(flds, 4);
        asserta(StrToUint(flds[0]) == FeatureIdx);
        const string &FeatureName = flds[1];
        uint AlphaSize = StrToUint(flds[2]);
        float Weight = (float) StrToFloat(flds[3]);
        
        m_FeatureNames.push_back(FeatureName);
        m_AlphaSizes.push_back(AlphaSize);
        m_Weights.push_back(Weight);
        
        vector<float> &LogProbs = m_LogProbsVec[FeatureIdx];
        GetNextFields(flds, AlphaSize + 1);
        asserta(flds[0] == "freqs");
        for (uint Letter = 0; Letter < AlphaSize; ++Letter)
        {
            // Probability is estimated as frequency
            float Prob = (float) StrToFloat(flds[Letter+1]);
            if (Prob < VERY_SMALL_FREQ)
                Prob = VERY_SMALL_FREQ;
            float LogProb = logf(Prob);
            LogProbs.push_back(LogProb);
        }
        
        vector<vector<float> > &LogProbMx = m_LogProbMxVec[FeatureIdx];
        LogProbMx.resize(AlphaSize);
        for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
            LogProbMx[Letter1].resize(AlphaSize);
        for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
        {
            GetNextFields(flds, Letter1 + 2);
            asserta(StrToUint(flds[0]) == Letter1);
            for (uint Letter2 = 0; Letter2 <= Letter1; ++Letter2)
            {
                // Probability is estimated as frequency
                float Prob = (float) StrToFloat(flds[Letter2+1]);
                if (Prob < VERY_SMALL_FREQ)
                    Prob = VERY_SMALL_FREQ;
                float LogProb = logf(Prob);
                LogProbMx[Letter1][Letter2] = LogProb;
                LogProbMx[Letter2][Letter1] = LogProb;
            }
        }
    }
//    ProfileCount = 2;
    m_Profiles.resize(ProfileCount);
    m_Seqs.resize(ProfileCount);

    for (uint ProfileIdx = 0; ProfileIdx < ProfileCount; ++ProfileIdx)
    {
        vector<vector<byte> > &Profile = m_Profiles[ProfileIdx];
        string &S = m_Seqs[ProfileIdx];
        GetNextFields(flds, 4);
        asserta(flds[0] == "chain");
        asserta(StrToUint(flds[1]) == ProfileIdx);
        const string &Label = flds[2];
        m_Labels.push_back(Label);
        const uint L = StrToUint(flds[3]);
        Profile.resize(L);
        for (uint Pos = 0; Pos < L; ++Pos)
        {
            GetNextFields(flds, 3);
            asserta(StrToUint(flds[0]) == ProfileIdx);
            asserta(StrToUint(flds[1]) == Pos);
            const string &Syms = flds[2];
            asserta(SIZE(Syms) == m_FeatureCount);
            for (uint FeatureIdx = 0; FeatureIdx < m_FeatureCount; ++FeatureIdx)
            {
                byte Sym = Syms[FeatureIdx];
                if (FeatureIdx == 0)
                {
                    uint Letter = g_CharToLetterAmino[Sym];
                    Profile[Pos].push_back(Letter);
                    S += Sym;
                }
                else
                {
                    uint Letter = uint(Sym - 'A');
                    asserta(Letter < 16);
                    Profile[Pos].push_back(Letter);
                }
            }
        }
    }
    m_Lines.clear();
//    m_FeatureCount = 1;
//    m_Weights[0] = 1.0;
/*
    bool Nucleo = (g_Alpha == ALPHA_Nucleo);
    HMMParams HP;
    HP.FromDefaults(Nucleo);
    HP.ToPairHMM();
    
    for (uint Letter1 = 0; Letter1 < 20; ++Letter1)
        for (uint Letter2 = 0; Letter2 < 20; ++Letter2)
        {
            byte c1 = g_LetterToCharAmino[Letter1];
            byte c2 = g_LetterToCharAmino[Letter2];
            m_LogProbMxVec[0][Letter1][Letter2] = PairHMM::m_MatchScore[c1][c2];
        }
    for (uint Letter1 = 0; Letter1 < 20; ++Letter1)
    {
        byte c1 = g_LetterToCharAmino[Letter1];
        m_LogProbsVec[0][Letter1] = PairHMM::m_InsScore[c1];
    }
*/
}
float Mega::GetInsScore(const vector<vector<byte> > &Profile, uint Pos) const
	{
	asserta(Pos < SIZE(Profile));
	const vector<byte> &ProfCol = Profile[Pos];
	float Score = 0;
	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		const vector<float> &LogProbs = m_LogProbsVec[i];
		byte Letter = ProfCol[i];
		Score += LogProbs[Letter]*m_Weights[i];
		}
	return Score;
	}

float Mega::GetMatchScore(
  const vector<vector<byte> > &ProfileX, uint PosX,
  const vector<vector<byte> > &ProfileY, uint PosY) const
	{
	asserta(PosX < SIZE(ProfileX));
	asserta(PosY < SIZE(ProfileY));
	const vector<byte> &ProfColX = ProfileX[PosX];
	const vector<byte> &ProfColY = ProfileY[PosY];
	float Score = 0;
	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		const vector<vector<float> > &SubstMx = m_LogProbMxVec[i];
		byte LetterX = ProfColX[i];
		byte LetterY = ProfColY[i];
		Score += SubstMx[LetterX][LetterY]*m_Weights[i];
		}
	return Score;
	}

void Mega::LogVec(const string &Name, const vector<float> &Vec) const
	{
	const uint N = SIZE(Vec);
	Log("\n%s/%u", Name.c_str(), N);
	for (uint i = 0; i < N; ++i)
		{
		if (i%10 == 0)
			Log("\n  ");
		else
			Log(" ");
		Log("[%2u]=%.2f", i, Vec[i]);
		}
	Log("\n");
	}

void Mega::LogMx(const string &Name,
  const vector<vector<float> > &Mx) const
	{
	const uint N = SIZE(Mx);
	Log("\n%s/%u\n", Name.c_str(), N);

	Log("     ");
	for (uint j = 0; j < N; ++j)
		Log(" %7u", j);
	Log("\n");
	for (uint i = 0; i < N; ++i)
		{
		Log("[%2u] ", i);
		const vector<float> &Row = Mx[i];
		asserta(SIZE(Row) == N);
		for (uint j = 0; j < N; ++j)
			Log(" %7.2f", Row[j]);
		Log("\n");
		}
	}

void Mega::LogFeatureParams(uint Idx) const
	{
	asserta(Idx < SIZE(m_FeatureNames));
	asserta(Idx < SIZE(m_LogProbMxVec));
	asserta(Idx < SIZE(m_LogProbsVec));
	const string &Name = m_FeatureNames[Idx];
	Log("\n");
	Log("Feature %s, weight %.3g\n",
	  Name.c_str(), m_Weights[Idx]);
	LogVec(Name, m_LogProbsVec[Idx]);
	LogMx(Name, m_LogProbMxVec[Idx]);
	}
