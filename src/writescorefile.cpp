#include "muscle.h"
#include "msa.h"
#include <errno.h>

extern float VTML_SP[32][32];
extern float NUC_SP[32][32];

static double GetColScore(const MSA &msa, unsigned uCol)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	unsigned uPairCount = 0;
	double dSum = 0.0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		if (msa.IsGap(uSeq1, uCol))
			continue;
		unsigned uLetter1 = msa.GetLetterEx(uSeq1, uCol);
		if (uLetter1 >= g_AlphaSize)
			continue;
		for (unsigned uSeq2 = uSeq1 + 1; uSeq2 < uSeqCount; ++uSeq2)
			{
			if (msa.IsGap(uSeq2, uCol))
				continue;
			unsigned uLetter2 = msa.GetLetterEx(uSeq2, uCol);
			if (uLetter2 >= g_AlphaSize)
				continue;
			double Score;
			switch (g_Alpha)
				{
			case ALPHA_Amino:
				Score = VTML_SP[uLetter1][uLetter2];
				break;
			case ALPHA_DNA:
			case ALPHA_RNA:
				Score = NUC_SP[uLetter1][uLetter2];
				break;
			default:
				Quit("GetColScore: invalid alpha=%d", g_Alpha);
				}
			dSum += Score;
			++uPairCount;
			}
		}
	if (0 == uPairCount)
		return 0;
	return dSum / uPairCount;
	}

void WriteScoreFile(const MSA &msa)
	{
	FILE *f = fopen(g_pstrScoreFileName, "w");
	if (0 == f)
		Quit("Cannot open score file '%s' errno=%d", g_pstrScoreFileName, errno);

	const unsigned uColCount = msa.GetColCount();
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uCol = 0; uCol < uColCount; ++uCol)
		{
		double Score = GetColScore(msa, uCol);
		fprintf(f, "%10.3f  ", Score);
		for (unsigned uSeq = 0; uSeq < uSeqCount; ++uSeq)
			{
			char c = msa.GetChar(uSeq, uCol);
			fprintf(f, "%c", c);
			}
		fprintf(f, "\n");
		}
	fclose(f);
	}
