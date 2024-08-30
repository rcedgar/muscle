#include "muscle.h"
#include "masmcol.h"
#include "masm.h"
#include "mega.h"
#include "sort.h"

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

void MASMCol::SetSortOrderVec()
	{
	const uint FeatureCount = SIZE(m_FreqsVec);
	m_SortOrderVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		vector<byte> &SortOrder = m_SortOrderVec[FeatureIdx];
		const vector<float> &Freqs = m_FreqsVec[FeatureIdx];
		vector<uint> Order(AlphaSize);
		QuickSortOrderDesc(Freqs.data(), AlphaSize, Order.data());
		for (uint i = 0; i < AlphaSize; ++i)
			{
			uint Letter = Order[i];
			asserta(Letter < AlphaSize);
			SortOrder.push_back(byte(Letter));
			}
		}
	}

void MASMCol::ToFile(FILE *f, uint ColIndex) const
	{
	if (f == 0)
		return;
	const uint FeatureCount = SIZE(m_FreqsVec);
	asserta(SIZE(m_ScoresVec) == FeatureCount);
	asserta(SIZE(m_SortOrderVec) == FeatureCount);
	fprintf(f, "col\t%u\t%u\n", ColIndex, FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
		const string &Name = Mega::GetFeatureName(FeatureIdx);
		fprintf(f, "feature\t%u\t%s\n", FeatureIdx, Name.c_str());

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

		const vector<byte> &SortOrder = m_SortOrderVec[FeatureIdx];
		fprintf(f, "sort");
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			fprintf(f, "\t%u", SortOrder[Letter]);
		fprintf(f, "\n");
		}
	}
