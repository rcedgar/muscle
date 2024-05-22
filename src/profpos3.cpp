#include "muscle.h"
#include "profpos3.h"
#include "alpha.h"
#include "m3alnparams.h"

uint ProfPos3::m_NewCount;
uint ProfPos3::m_DeleteCount;

void SortCounts(const float *Counts, byte *Order)
	{
	static byte InitialSortOrder[20] =
		{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
		};
	memcpy(Order, InitialSortOrder, g_AlphaSize*sizeof(byte));

	bool bAny = true;
	while (bAny)
		{
		bAny = false;
		for (unsigned n = 0; n < g_AlphaSize - 1; ++n)
			{
			unsigned i1 = Order[n];
			unsigned i2 = Order[n+1];
			if (Counts[i1] < Counts[i2])
				{
				Order[n+1] = i1;
				Order[n] = i2;
				bAny = true;
				}
			}
		}
	}

// Pseudo-column before first actual column
void ProfPos3::SetStartDimers()
	{
	m_LL = 1.0f;
	m_LG = 0.0f;
	m_GL = 0.0f;
	m_GG = 0.0f;
	}

void ProfPos3::SetOcc()
	{
	m_fOcc = m_LL + m_GL;
	}

// Set alignment scores from frequencies
void ProfPos3::SetAAScores(const Mx2020 &SubstMx_Letter)
	{
	SortCounts(m_Freqs, m_SortOrder);
	for (unsigned i = 0; i < g_AlphaSize; ++i)
		{
		float Sum = 0;
		for (unsigned j = 0; j < g_AlphaSize; ++j)
			{
			float Freq = m_Freqs[j];
			float Score = SubstMx_Letter[i][j];
			Sum += Freq*Score;
			}
		m_AAScores[i] = Sum;
		}
	}

// SeqWeights must sum to 1 (for ProfPos3.m_GapOpen/Close)
void ProfPos3::SetFreqs(const MultiSequence &MSA, uint ColIndex,
  const vector<float> &SeqWeights)
	{
	const uint SeqCount = MSA.GetSeqCount();
#if DEBUG
	{
	assert(SIZE(SeqWeights) == SeqCount);
	float SumWeights = 0;
	for (uint i = 0; i < SeqCount; ++i)
		SumWeights += SeqWeights[i];
	assert(SumWeights > 0.9 && SumWeights < 1.1);
	}
#endif
	const uint ColCount = MSA.GetColCount();
	const uint LastCol = ColCount - 1;

	m_AllGaps = true;
	memset_zero(m_Freqs, 20);
	m_LL = 0;
	m_LG = 0;
	m_GL = 0;
	m_GG = 0;
	m_fOcc = 0;

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		float w = SeqWeights[SeqIndex];
		const byte *ByteSeq = MSA.GetBytePtr(SeqIndex);
		byte c = ByteSeq[ColIndex];
		bool LetterHere = !isgap(c);
		bool LetterPrev = (ColIndex == 0 || !isgap(ByteSeq[ColIndex-1]));
		bool LetterNext = (ColIndex == LastCol || !isgap(ByteSeq[ColIndex+1]));
		if (LetterHere)
			{
			m_AllGaps = false;
			m_fOcc += w;
			uint Letter = g_CharToLetter[c];
			if (Letter < g_AlphaSize)
				m_Freqs[Letter] += w;

			if (LetterPrev)
				m_LL += w;
			else
				m_GL += w;
			}
		else
			{
			if (LetterPrev)
				m_LG += w;
			else
				m_GG += w;
			}
		}
	//NormalizeAAFreqsIfRequired();
	}

void ProfPos3::SetFreqs2(uint SeqCount,
  uint LLCount, uint LGCount, uint GLCount, uint GGCount,
  const vector<uint> &LetterCounts)
	{
	asserta(SIZE(LetterCounts) == g_AlphaSize);
#if DEBUG
	{
	assert(LLCount + LGCount + GLCount + GGCount == SeqCount);
	uint LetterCount = LLCount + GLCount;
	uint LetterSum = 0;
	for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
		LetterSum += LetterCounts[Letter];
	assert(LetterSum <= SeqCount);
	assert(LetterSum == LetterCount);
	}
#endif
	float N = float(SeqCount);
	m_LL = float(LLCount)/N;
	m_LG = float(LGCount)/N;
	m_GL = float(GLCount)/N;
	m_GG = float(GGCount)/N;
	m_fOcc = m_LL + m_GL;

	for (uint Letter = 0; Letter < 20; ++Letter)
		m_Freqs[Letter] = 0;

	float SumFreq = 0.0f;
	for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
		{
		uint n = LetterCounts[Letter];
		float Freq = float(n)/N;
		m_Freqs[Letter] = Freq;
		SumFreq += Freq;
		}
	asserta(feq(SumFreq, m_fOcc));

	//NormalizeAAFreqsIfRequired();
	}

//void ProfPos3::NormalizeAAFreqsIfRequired()
//	{
//	//if (!g_NormalizeAAFreqs)
//	//	return;
//
//	float Sum = 0;
//	for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
//		Sum += m_Freqs[Letter];
//	if (Sum == 0)
//		return;
//
//	for (uint Letter = 0; Letter < g_AlphaSize; ++Letter)
//		m_Freqs[Letter] /= Sum;
//	}

void ProfPos3::LogMe() const
	{
#define w(x)	Log(" %s=%.3f", #x, m_##x)
	w(LL);
	w(LG);
	w(GL);
	w(GG);
	w(fOcc);
	//w(GapOpenScore);
	//w(GapCloseScore);
#undef w
	Log("\n");

	Log(" Freqs: ");
	bool ZeroFound = false;
	for (uint i = 0; i < 20; ++i)
		{
		uint Letter = m_SortOrder[i];
		float Freq = m_Freqs[Letter];
		char c = g_LetterToChar[Letter];
		if (Freq == 0)
			ZeroFound = true;
		if (ZeroFound)
			{
			asserta(Freq == 0);
			continue;
			}
		Log(" %c=%.3g", c, Freq);
		}
	Log("\n");

	Log("Scores: ");
	for (uint i = 0; i < 20; ++i)
		{
		if (i == 9)
			Log("\n   ");
		uint Letter = m_SortOrder[i];
		float Score = m_AAScores[Letter];
		char c = g_LetterToChar[Letter];
		Log(" %c=%.3g", c, Score);
		}
	Log("\n");
	}

void ProfPos3::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;

#define w(x)	fprintf(f, "\t%.5g", m_##x)

	w(LL);
	w(LG);
	w(GL);
	w(GG);
	w(fOcc);
	w(GapOpenScore);
	w(GapCloseScore);

	for (int i = 0; i < 20; ++i)
		{
		w(Freqs[i]);
		w(AAScores[i]);
		}

#undef w

	fprintf(f, "\n");
	}
