#include "muscle.h"
#include "profile3.h"
#include "m3alnparams.h"

uint Profile3::m_NewCount;
uint Profile3::m_DeleteCount;

void Profile3::Clear()
	{
	if (m_PPs.empty())
		return;
	size_t N = m_PPs.size();
	for (size_t i = 0; i < N; ++i)
		delete m_PPs[i];
	m_PPs.clear();
	}

const ProfPos3 &Profile3::GetPP(uint ColIndex) const
	{
	assert(ColIndex < SIZE(m_PPs));
	return *m_PPs[ColIndex];
	}

void Profile3::SetGapOpenScore(float GapOpen, uint ColIndex)
	{
	assert(ColIndex < SIZE(m_PPs));
	ProfPos3 *PP = m_PPs[ColIndex];
	if (ColIndex == 0)
		PP->m_GapOpenScore = PP->m_fOcc*GapOpen/2;
	else
		{
		float GapOpenFreq = PP->m_LG;
		PP->m_GapOpenScore = GapOpen*(1.0f - GapOpenFreq)/2;
		}
	}

void Profile3::SetGapCloseScore(float GapOpen, uint ColIndex)
	{
	const uint ColCount = GetColCount();
	assert(ColIndex < ColCount);
	ProfPos3 *PP = m_PPs[ColIndex];
	if (ColIndex + 1 == ColCount)
		PP->m_GapCloseScore = GapOpen*PP->m_fOcc/2;
	else
		{
		const ProfPos3 *NextPP = m_PPs[ColIndex+1];
		float GapCloseFreq = NextPP->m_GL;
		PP->m_GapCloseScore = GapOpen*(1.0f - GapCloseFreq)/2;
		}
	}

void Profile3::FromSeq(const Sequence &Seq,
  const Mx2020 &SubstMx_Letter, float GapOpen)
	{
	MultiSequence MSA;
	MSA.AddSequence(&Seq, false);
	vector<float> SeqWeights;
	SeqWeights.push_back(1.0f);
	FromMSA(MSA, SubstMx_Letter, GapOpen, SeqWeights);
	}

// SeqWeights must sum to 1 (for ProfPos3.m_GapOpen/Close)
void Profile3::FromMSA(const MultiSequence &MSA,
  const Mx2020 &SubstMx_Letter, float GapOpen,
  const vector<float> &SeqWeights)
	{
	const uint ColCount = MSA.GetColCount();
	const uint SeqCount = MSA.GetSeqCount();
	asserta(SIZE(SeqWeights) == SeqCount);
#if 1 // DEBUG
	{
	float SumWeights = 0;
	for (uint i = 0; i < SeqCount; ++i)
		SumWeights += SeqWeights[i];
	asserta(SumWeights > 0.9 && SumWeights < 1.1);
	}
#endif
	m_PPs.clear();
	m_PPs.reserve(ColCount);
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		ProfPos3 *PP = new ProfPos3;
		PP->SetFreqs(MSA, ColIndex, SeqWeights);
		PP->SetAAScores(SubstMx_Letter);
		m_PPs.push_back(PP);
		}
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		ProfPos3 *PP = m_PPs[ColIndex];
		SetGapOpenScore(GapOpen, ColIndex);
		SetGapCloseScore(GapOpen, ColIndex);
		}
	}

static void LogF(float f)
	{
	if (f > -0.00001 && f < 0.00001)
		Log("       ");
	else
		Log("  %5.3f", f);
	}

void Profile3::LogMe(const MultiSequence *MSA) const
	{
	Log("  Pos  Occ     LL     LG     GL     GG       Open    Close\n");
	Log("  ---  ---     --     --     --     --     ------  -------\n");
	const uint ColCount = GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const ProfPos3 &PP = *m_PPs[ColIndex];
		Log("%5u", ColIndex);
		LogF(PP.m_fOcc);
		LogF(PP.m_LL);
		LogF(PP.m_LG);
		LogF(PP.m_GL);
		LogF(PP.m_GG);
		Log("  %7.3f", -PP.m_GapOpenScore);
		Log("  %7.3f", -PP.m_GapCloseScore);
		if (MSA != 0)
			{
			const uint uSeqCount = MSA->GetSeqCount();
			Log("  ");
			for (uint uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", MSA->GetChar(uSeqIndex, ColIndex));
			}
		Log("\n");
		}

	Log("\n");
	Log("  Pos");
	for (uint n = 0; n < g_AlphaSize; ++n)
		Log("     %c", g_LetterToChar[n]);
	Log("\n");
	Log("  ---");
	for (uint n = 0; n < g_AlphaSize; ++n)
		Log(" -----");
	Log("\n");

	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const ProfPos3 &PP = *m_PPs[ColIndex];
		Log("%5u", ColIndex);

		float SumFreqs = 0;
		for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
			{
			float f = PP.m_Freqs[Letter];
			SumFreqs += f;
			if (f == 0.0)
				Log("      ");
			else
				Log(" %5.3f", f);
			}
		Log("  (Sum=%4.2f)", SumFreqs);
		if (MSA != 0)
			{
			const uint uSeqCount = MSA->GetSeqCount();
			Log("  ");
			for (uint uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", MSA->GetChar(uSeqIndex, ColIndex));
			}
		for (uint k = 0; k < g_AlphaSize; ++k)
			{
			uint i = PP.m_SortOrder[k];
			float Freq = PP.m_Freqs[i];
			if (Freq > 0)
				Log(" %c(%.3f)", g_LetterToChar[i], Freq);
			}

		Log("\n");
		}

	Log("\n");
	Log("Scores\n");
	Log("  Col");
	for (uint n = 0; n < g_AlphaSize; ++n)
		Log("     %c", g_LetterToChar[n]);
	Log("\n");
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const ProfPos3 &PP = *m_PPs[ColIndex];
		Log("%5u", ColIndex);

		float SumFreqs = 0;
		for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
			{
			float f = PP.m_AAScores[Letter];
			SumFreqs += f;
			Log(" %5.2f", f);
			}
		Log("\n");
		}
	}

void Profile3::SetAAScores(const Mx2020 &SubstMx_Letter)
	{
	const uint ColCount = SIZE(m_PPs);
	for (uint Col = 0; Col < ColCount; ++Col)
		m_PPs[Col]->SetAAScores(SubstMx_Letter);
	}

void Profile3::SetScores(const Mx2020 &SubstMx_Letter, float GapOpen)
	{
	SetAAScores(SubstMx_Letter);
	SetGapScores(GapOpen);
	}

void Profile3::SetGapScores(float GapOpen)
	{
	const uint ColCount = GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		SetGapOpenScore(GapOpen, ColIndex);
		SetGapCloseScore(GapOpen, ColIndex);
		}
	}

void Profile3::Validate() const
	{
	const uint ColCount = GetColCount();
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const ProfPos3 &PP = *m_PPs[ColIndex];
		if (!feq(PP.m_fOcc, PP.m_LL + PP.m_GL))
			Die("Col %u, fOcc != LL + GL", ColIndex);

		float s1 = PP.m_LL + PP.m_LG + PP.m_GL + PP.m_GG;
		asserta(feq(s1, 1.0));

		if (ColIndex > 0)
			{
			const ProfPos3 &PPPrev = *m_PPs[ColIndex-1];
			float s2 = PPPrev.m_LL + PPPrev.m_GL;
			float s3 = PP.m_LL + PP.m_LG;
			if (!feq(s2, s3))
				Die("Col %u, LL + LG != Prev.LL + Prev.GL", ColIndex);
			}
		if (ColIndex + 1 < ColCount)
			{
			const ProfPos3 &PPNext = *m_PPs[ColIndex+1];
			float s4 = PP.m_LL + PP.m_GL;
			float s5 = PPNext.m_LL + PPNext.m_LG;
			if (!feq(s4, s5))
				Die("Col %u, LL + GL != Next.LL + Next.LG", ColIndex);
			}
		}
	}

void Profile3::ToTsv(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToTsv(f);
	CloseStdioFile(f);
	}

void Profile3::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	const uint ColCount = GetColCount();
	fprintf(f, "%u\n", ColCount);
	for (uint i = 0; i < ColCount; ++i)
		{
		fprintf(f, "%u", i);
		m_PPs[i]->ToTsv(f);
		}
	}

float Profile3::GetSelfScore() const
	{
	float Sum = 0;
	const uint ColCount = GetColCount();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		const ProfPos3 &PP = GetPP(Col);
		float Score = ScoreProfPos2(PP, PP);
		Sum += Score;
		}
	return Sum;
	}

uint Profile3::LogDiffs(const Profile3 &Prof2) const
	{
	const uint ColCount = GetColCount();
	const uint ColCount2 = Prof2.GetColCount();
	if (ColCount2 != ColCount)
		{
		Log("Lengths differ %u, %u\n", ColCount, ColCount2);
		return 1;
		}

	uint DiffCount = 0;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const ProfPos3 &PP = GetPP(ColIndex);
		const ProfPos3 &PP2 = Prof2.GetPP(ColIndex);

#define w(x)	if (!feq(PP.m_##x, PP2.m_##x)) \
					{ \
					++DiffCount; \
					Log("Col %u %s: %.4g %.4g\n", ColIndex, #x, PP.m_##x, PP2.m_##x); \
					}

			w(LL);
			w(LG);
			w(GL);
			w(GG);
			w(fOcc);
			w(GapOpenScore);
			w(GapCloseScore);

#undef w

			for (int i = 0; i < 20; ++i)
				{
				float f = PP.m_Freqs[i];
				float f2 = PP2.m_Freqs[i];
				if (!feq(f, f2))
					{
					++DiffCount;
					Log("Col %u Freqs[%u=%c] %.4g %.4g\n",
					  ColIndex, i, g_LetterToChar[i], f, f2);
					}

				f = PP.m_AAScores[i];
				f2 = PP2.m_AAScores[i];
				if (!feq(f, f2))
					{
					++DiffCount;
					Log("Col %u AAScores[%u=%c] %.4g %.4g\n",
					  ColIndex, i, g_LetterToChar[i], f, f2);
					}
				}
		}
	return DiffCount;
	}

void cmd_build_prof3()
	{
	const string &InputFileName = opt(build_prof3);

	MultiSequence MSA;
	MSA.FromFASTA(InputFileName);
	bool IsNucleo = MSA.GuessIsNucleo();
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);
	M3AlnParams AP;
	AP.SetFromCmdLine(IsNucleo);
	const Mx2020 &SubstMx_Letter =
	  AP.m_SubstMx_Letter;
	float GapOpen = AP.m_GapOpen;

	const uint SeqCount = MSA.GetSeqCount();
	float w = 1.0f/SeqCount;
	vector<float> SeqWeights(SeqCount, w);

	Profile3 Prof;
	Prof.FromMSA(MSA, SubstMx_Letter, GapOpen, SeqWeights);
	Prof.LogMe(&MSA);
	Prof.Validate();
	Prof.ToTsv(opt(output));
	}
