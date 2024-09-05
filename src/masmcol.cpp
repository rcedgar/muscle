#include "muscle.h"
#include "masmcol.h"
#include "masm.h"
#include "mega.h"
#include "sort.h"

char MASMCol::GetConsensusAAChar() const
	{
	uint FI = m_MASM->m_AAFeatureIdx;
	asserta(FI < SIZE(m_FreqsVec));
	const vector<float> &Freqs = m_FreqsVec[FI];
	asserta(SIZE(Freqs) == 20);
	float MaxFreq = 0;
	float SumFreq = 0;
	uint MaxLetter = 0;
	for (uint Letter = 0; Letter < 20; ++Letter)
		{
		float Freq = Freqs[Letter];
		if (Freq > MaxFreq)
			{
			MaxFreq = Freq;
			MaxLetter = Letter;
			}
		SumFreq += Freq;
		}
	if (SumFreq < 0.5f)
		return '-';
	char c = g_LetterToCharAmino[MaxLetter];
	if (MaxFreq < 0.5)
		return tolower(c);
	return c;
	}

void MASMCol::SetScoreVec()
	{
	const uint FeatureCount = SIZE(m_FreqsVec);
	m_ScoresVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);

	// Weights are already applied to ScoreMx
		const vector<vector<float> > &ScoreMx = Mega::m_LogOddsMxVec[FeatureIdx];
		asserta(SIZE(ScoreMx) == AlphaSize);
		const vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		asserta(SIZE(Freqs) == AlphaSize);
		vector<float> &Scores = m_ScoresVec[FeatureIdx];
		for (byte Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float Total = 0;
			for (byte Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
				{
				float Freq2 = Freqs[Letter2];
				Total += Freq2*ScoreMx[Letter][Letter2];
				}

		// Weights are already applied to ScoreMx
			Scores.push_back(Total);
			}
		}
	}

const vector<float> &MASMCol::GetAAScores() const
	{
	uint AAFeatureIdx = m_MASM->m_AAFeatureIdx;
	asserta(AAFeatureIdx < SIZE(m_ScoresVec));
	const vector<float> &Scores = m_ScoresVec[AAFeatureIdx];
	return Scores;
	}

void MASMCol::LogMe() const
	{
	Log("  MSAMCol[%u]", m_ColIndex);
	const uint FeatureCount = SIZE(m_FreqsVec);
	if (FeatureCount == 1)
		{
		const uint FeatureIdx = 0;
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		const string &Name = Mega::GetFeatureName(FeatureIdx);
		const vector<float> &Scores = m_ScoresVec[FeatureIdx];
		Log("  | %s |", Name.c_str());
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Log(" %c=%.3g", g_LetterToCharAmino[Letter], Scores[Letter]);
		Log("\n");
		return;
		}

	Log("\n");
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		const string &Name = Mega::GetFeatureName(FeatureIdx);
		const vector<float> &Scores = m_ScoresVec[FeatureIdx];
		Log("  | %8.8s |", Name.c_str());
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Log(" %c=%.3g", g_LetterToCharAmino[Letter], Scores[Letter]);
		Log("\n");
		}
	}

void MASMCol::ToFile(FILE *f) const
	{
	if (f == 0)
		return;
	const uint FeatureCount = SIZE(m_FreqsVec);
	asserta(SIZE(m_ScoresVec) == FeatureCount);
	//asserta(SIZE(m_SortOrderVec) == FeatureCount);
	fprintf(f, "col\t%u\n", m_ColIndex);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		const string &Name = Mega::GetFeatureName(FeatureIdx);
		fprintf(f, "colfeature\t%u\n", FeatureIdx);

		const vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		fprintf(f, "freqs");
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			fprintf(f, "\t%.3g", Freqs[Letter]);
		fprintf(f, "\n");

		const vector<float> &Scores = m_ScoresVec[FeatureIdx];
		fprintf(f, "scores");
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			fprintf(f, "\t%.3g", Scores[Letter]);
		fprintf(f, "\n");

		//const vector<byte> &SortOrder = m_SortOrderVec[FeatureIdx];
		//fprintf(f, "sort");
		//for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		//	fprintf(f, "\t%u", SortOrder[Letter]);
		//fprintf(f, "\n");
		}
	}

void MASMCol::FromFile(FILE *f)
	{
	asserta(f != 0);
	asserta(m_MASM != 0);
	const uint FeatureCount = m_MASM->m_FeatureCount;
	vector<string> Fields;
	ReadTabbedLine(f, Fields, 2);
	asserta(Fields[0] == "col");
	m_ColIndex = StrToUint(Fields[1]);
	m_FreqsVec.clear();
	m_ScoresVec.clear();
	m_FreqsVec.resize(FeatureCount);
	m_ScoresVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		ReadTabbedLine(f, Fields, 2);
		asserta(Fields[0] == "colfeature");
		asserta(StrToUint(Fields[1]) == FeatureIdx);
		uint AlphaSize = m_MASM->m_AlphaSizes[FeatureIdx];

		vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		vector<float> &Scores = m_ScoresVec[FeatureIdx];
		ReadTabbedLine(f, Fields, AlphaSize+1);
		asserta(Fields[0] == "freqs");
		float SumFreqs = 0;
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			{
			float Freq = (float) StrToFloat(Fields[Letter+1]);
			Freqs.push_back(Freq);
			SumFreqs += Freq;
			}
		asserta(SumFreqs < 1.001);

		ReadTabbedLine(f, Fields, AlphaSize+1);
		asserta(Fields[0] == "scores");
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Scores.push_back((float) StrToFloat(Fields[Letter+1]));
		}
	}

float MASMCol::GetMatchScore_MegaProfilePos(const vector<byte> &ProfPos) const
	{
	float Total = 0;
	const uint FeatureCount = SIZE(ProfPos);
	assert(FeatureCount == m_MASM->m_FeatureCount);
	assert(SIZE(m_ScoresVec) == FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		byte Letter = ProfPos[FeatureIdx];
		Total += m_ScoresVec[FeatureIdx][Letter];
		}
	return Total;
	}
