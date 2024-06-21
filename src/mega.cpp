#include "myutils.h"
#include "mega.h"
#include "alpha.h"

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

void Mega::CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
  vector<float> &Freqs) const
	{
	AssertSymmetrical(FreqsMx);
	const uint N = SIZE(FreqsMx);
	}

void Mega::CalcLogProbsMx(const vector<vector<float > > &FreqsMx,
  vector<vector<float > > &LogProbMx)
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
	uint m_FeatureCount = StrToUint(flds[1]);
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

		vector<float> &Freqs = m_LogProbsVec[FeatureIdx];
		GetNextFields(flds, AlphaSize + 1);
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Freqs.push_back((float) StrToFloat(flds[Letter+1]));

		vector<vector<float> > &LogProbMx = m_LogProbMxVec[FeatureIdx];
		LogProbMx.resize(AlphaSize);
		for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
			LogProbMx[Letter1].resize(AlphaSize, FLT_MAX);

		for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
			{
			GetNextFields(flds, Letter1 + 2);
			asserta(StrToUint(flds[0]) == Letter1);
			for (uint Letter2 = 0; Letter2 <= Letter1; ++Letter2)
				{
				float Score = (float) log(StrToFloat(flds[Letter2+1]));
				LogProbMx[Letter1][Letter2] = Score;
				LogProbMx[Letter2][Letter1] = Score;
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
		Score += LogProbs[i]*m_Weights[i];
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

void Mega::LogMx(const string &Name,
  const vector<vector<float> > &Mx) const
	{
	const uint N = SIZE(Mx);
	Log("\n%s(%u)\n", Name.c_str(), N);

	Log("     ");
	for (uint j = 0; j < N; ++j)
		Log(" %10u", j);
	Log("\n");
	for (uint i = 0; i < N; ++i)
		{
		Log("[%2u] ", i);
		const vector<float> &Row = Mx[i];
		asserta(SIZE(Row) == N);
		for (uint j = 0; j < N; ++j)
			Log(" %10.3g", Row[j]);
		Log("\n");
		}
	}

void Mega::LogFeatureParams(uint Idx) const
	{
	asserta(Idx < SIZE(m_FeatureNames));
	asserta(Idx < SIZE(m_FreqsMxVec));
	asserta(Idx < SIZE(m_LogProbMxVec));
	asserta(Idx < SIZE(m_LogProbsVec));
	Log("\n");
	Log("Feature %s, weight %.3g\n",
	  m_FeatureNames[Idx].c_str(), m_Weights[Idx]);
	LogMx("Freqs", m_FreqsMxVec[Idx]);
	LogMx("LogProbs", m_LogProbMxVec[Idx]);
	}
