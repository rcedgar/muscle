#include "muscle.h"
#include "msa.h"
#include "distfunc.h"
#include "seq.h"
#include "msa.h"

extern int BLOSUM62[20][20];
extern double BLOSUM62_Expected;

double ScoreDistSigma(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2,
  unsigned *ptrLength);

double GetMyDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2)
	{
	if (g_Alpha != ALPHA_Amino)
		Quit("GetMyDist() is only for amino acid sequences");

	const double GAP_OPEN = -5;
	const double GAP_EXT = -1;
	unsigned Length = 0;
	double TotalScore = 0;
	bool InGap = false;
	const unsigned ColCount = msa.GetColCount();
	const char *seq1 = msa.GetSeqBuffer(SeqIndex1);
	const char *seq2 = msa.GetSeqBuffer(SeqIndex2);
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c1 = seq1[ColIndex];
		char c2 = seq2[ColIndex];
		unsigned Letter1 = g_CharToLetterEx[c1];
		unsigned Letter2 = g_CharToLetterEx[c2];
		if (Letter1 == AX_GAP && Letter2 == AX_GAP)
			continue;

		++Length;
		if (Letter1 == AX_GAP || Letter2 == AX_GAP)
			{
			if (InGap)
				TotalScore += GAP_EXT;
			else
				{
				InGap = true;
				TotalScore += GAP_OPEN;
				}
			}
		else
			InGap = false;

		if (Letter1 >= 20 || Letter2 >= 20)
			continue;
		TotalScore += BLOSUM62[Letter1][Letter2];
		}

	if (Length == 0)
		return GAP_OPEN;

	double AvgScore = TotalScore/Length;
	return AvgScore;
	}
