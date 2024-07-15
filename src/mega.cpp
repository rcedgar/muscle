#include "myutils.h"
#include "mega.h"
#include "alpha.h"
#include "pairhmm.h"

static const float VERY_SMALL_FREQ = 1e-6f;

string Mega::m_FileName;
vector<string> Mega::m_Lines;
vector<string> Mega::m_FeatureNames;
vector<float> Mega::m_Weights;
vector<uint> Mega::m_AlphaSizes;
vector<string> Mega::m_Labels;
vector<vector<vector<byte> > > Mega::m_Profiles;
vector<string> Mega::m_Seqs;
vector<vector<float> > Mega::m_LogProbsVec;
vector<vector<vector<float> > > Mega::m_LogProbMxVec;
vector<string> m_Labels;
uint Mega::m_NextLineNr;
uint Mega::m_FeatureCount;
unordered_map<string, uint> Mega::m_LabelToIdx;

uint Mega::GetGSIByLabel(const string &Label)
	{
	unordered_map<string, uint>::const_iterator iter = m_LabelToIdx.find(Label);
	if (iter == m_LabelToIdx.end())
		Die("Mega::GetGSIByLabel(%s)", Label.c_str());
	uint GSI = iter->second;
	return GSI;
	}

const string &Mega::GetLabelByGSI(uint GSI)
	{
	asserta(GSI < SIZE(m_Labels));
	return m_Labels[GSI];
	}

const vector<vector<byte> > *Mega::GetProfileByGSI(uint GSI)
	{
	asserta(GSI < SIZE(m_Profiles));
	return &m_Profiles[GSI];
	}

const vector<vector<byte> > *Mega::GetProfileByLabel(const string &Label)
	{
	unordered_map<string, uint>::const_iterator iter = m_LabelToIdx.find(Label);
	if (iter == m_LabelToIdx.end())
		Die("Mega::GetProfileByLabel(%s)", Label.c_str());
	uint Idx = iter->second;
	asserta(Idx < SIZE(m_Profiles));
	return &m_Profiles[Idx];
	}

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
	asserta(m_FeatureNames.empty());
	asserta(m_FeatureCount == 0);
	asserta(m_Profiles.empty());

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
		if (m_LabelToIdx.find(Label) != m_LabelToIdx.end())
			Die("Duplicate label in mega file >%s", Label.c_str());
		m_LabelToIdx[Label] = ProfileIdx;
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
