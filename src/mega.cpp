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
	uint FeatureCount = StrToUint(flds[1]);
	uint ProfileCount = StrToUint(flds[2]);
	m_FreqsVec.resize(FeatureCount);
	m_SubstMxVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		GetNextFields(flds, 4);
		asserta(StrToUint(flds[0]) == FeatureIdx);
		const string &FeatureName = flds[1];
		uint AlphaSize = StrToUint(flds[2]);
		float Weight = (float) StrToFloat(flds[3]);

		m_FeatureNames.push_back(FeatureName);
		m_AlphaSizes.push_back(AlphaSize);
		m_Weights.push_back(Weight);

		vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		GetNextFields(flds, AlphaSize + 1);
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Freqs.push_back((float) StrToFloat(flds[Letter+1]));

		vector<vector<float> > &SubstMx = m_SubstMxVec[FeatureIdx];
		SubstMx.resize(AlphaSize);
		for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
			SubstMx[Letter1].resize(AlphaSize);

		for (uint Letter1 = 0; Letter1 < AlphaSize; ++Letter1)
			{
			GetNextFields(flds, Letter1 + 2);
			asserta(StrToUint(flds[0]) == Letter1);
			for (uint Letter2 = 0; Letter2 <= Letter1; ++Letter2)
				{
				float Score = (float) StrToFloat(flds[Letter2+1]);
				SubstMx[Letter1][Letter2] = Score;
				SubstMx[Letter2][Letter1] = Score;
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
			asserta(SIZE(Syms) == FeatureCount);
			for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
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
