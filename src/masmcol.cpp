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
			Scores.push_back(Total);
			}
		}
	}

//void MASMCol::SetSortOrderVec()
//	{
//	const uint FeatureCount = SIZE(m_FreqsVec);
//	m_SortOrderVec.resize(FeatureCount);
//	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
//		{
//		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
//		vector<byte> &SortOrder = m_SortOrderVec[FeatureIdx];
//		const vector<float> &Freqs = m_FreqsVec[FeatureIdx];
//		vector<uint> Order(AlphaSize);
//		QuickSortOrderDesc(Freqs.data(), AlphaSize, Order.data());
//		for (uint i = 0; i < AlphaSize; ++i)
//			{
//			uint Letter = Order[i];
//			asserta(Letter < AlphaSize);
//			SortOrder.push_back(byte(Letter));
//			}
//		}
//	}

const vector<float> &MASMCol::GetAAScores() const
	{
	uint AAFeatureIdx = m_MASM->m_AAFeatureIdx;
	asserta(AAFeatureIdx < SIZE(m_ScoresVec));
	const vector<float> &Scores = m_ScoresVec[AAFeatureIdx];
	return Scores;
	}

void MASMCol::ToFile(FILE *f, uint ColIndex) const
	{
	if (f == 0)
		return;
	const uint FeatureCount = SIZE(m_FreqsVec);
	asserta(SIZE(m_ScoresVec) == FeatureCount);
	//asserta(SIZE(m_SortOrderVec) == FeatureCount);
	fprintf(f, "col\t%u\n", ColIndex);
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

void MASMCol::FromFile(FILE *f, uint ColIndex)
	{
	asserta(f != 0);
	asserta(m_MASM != 0);
	const uint FeatureCount = m_MASM->m_FeatureCount;
	//asserta(SIZE(m_SortOrderVec) == FeatureCount);
	//fprintf(f, "col\t%u\t%u\n", ColIndex, FeatureCount);
	vector<string> Fields;
	ReadTabbedLine(f, Fields, 2);
	asserta(Fields[0] == "col");
	asserta(StrToUint(Fields[1]) == ColIndex);
	m_FreqsVec.clear();
	m_ScoresVec.clear();
	m_FreqsVec.resize(FeatureCount);
	m_ScoresVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		//uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		//const string &Name = Mega::GetFeatureName(FeatureIdx);
		//fprintf(f, "feature\t%u\t%s\n", FeatureIdx, Name.c_str());
		ReadTabbedLine(f, Fields, 2);
		asserta(Fields[0] == "colfeature");
		asserta(StrToUint(Fields[1]) == FeatureIdx);
		uint AlphaSize = m_MASM->m_AlphaSizes[FeatureIdx];

		vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		vector<float> &Scores = m_ScoresVec[FeatureIdx];
		//fprintf(f, "freqs");
		//for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		//	fprintf(f, "\t%.3g", Freqs[Letter]);
		//fprintf(f, "\n");
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

		//const vector<float> &Scores = m_ScoresVec[FeatureIdx];
		//fprintf(f, "scores");
		//for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		//	fprintf(f, "\t%.3g", Scores[Letter]);
		//fprintf(f, "\n");
		ReadTabbedLine(f, Fields, AlphaSize+1);
		asserta(Fields[0] == "scores");
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			Scores.push_back((float) StrToFloat(Fields[Letter+1]));

		//const vector<byte> &SortOrder = m_SortOrderVec[FeatureIdx];
		//fprintf(f, "sort");
		//for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		//	fprintf(f, "\t%u", SortOrder[Letter]);
		//fprintf(f, "\n");
		}
	}
