#include "muscle.h"
#include "profile.h"

#define TRACE	0

enum
	{
	LL = 0,
	LG = 1,
	GL = 2,
	GG = 3,
	};

static const char *GapTypeToStr(int GapType)
	{
	switch (GapType)
		{
	case LL: return "LL";
	case LG: return "LG";
	case GL: return "GL";
	case GG: return "GG";
		}
	Quit("Invalid gap type");
	return "?";
	}

static SCORE GapScoreMatrix[4][4];

static void InitGapScoreMatrix()
	{
	const SCORE t = (SCORE) 0.2;

	GapScoreMatrix[LL][LL] = 0;
	GapScoreMatrix[LL][LG] = g_scoreGapOpen;
	GapScoreMatrix[LL][GL] = 0;
	GapScoreMatrix[LL][GG] = 0;

	GapScoreMatrix[LG][LL] = g_scoreGapOpen;
	GapScoreMatrix[LG][LG] = 0;
	GapScoreMatrix[LG][GL] = g_scoreGapOpen;
	GapScoreMatrix[LG][GG] = t*g_scoreGapOpen;	// approximation!

	GapScoreMatrix[GL][LL] = 0;
	GapScoreMatrix[GL][LG] = g_scoreGapOpen;
	GapScoreMatrix[GL][GL] = 0;
	GapScoreMatrix[GL][GG] = 0;

	GapScoreMatrix[GG][LL] = 0;
	GapScoreMatrix[GG][LG] = t*g_scoreGapOpen;	// approximation!
	GapScoreMatrix[GG][GL] = 0;
	GapScoreMatrix[GG][GG] = 0;

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < i; ++j)
			if (GapScoreMatrix[i][j] != GapScoreMatrix[j][i])
				Quit("GapScoreMatrix not symmetrical");
	}

static SCORE SPColBrute(const MSA &msa, unsigned uColIndex)
	{
	SCORE Sum = 0;
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		unsigned uLetter1 = msa.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter1 >= 20)
			continue;
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			unsigned uLetter2 = msa.GetLetterEx(uSeqIndex2, uColIndex);
			if (uLetter2 >= 20)
				continue;
			SCORE t = w1*w2*(*g_ptrScoreMatrix)[uLetter1][uLetter2];
#if	TRACE
			Log("Check %c %c w1=%.3g w2=%.3g Mx=%.3g t=%.3g\n",
			  LetterToCharAmino(uLetter1),
			  LetterToCharAmino(uLetter2),
			  w1,
			  w2,
			  (*g_ptrScoreMatrix)[uLetter1][uLetter2],
			  t);
#endif
			Sum += t;
			}
		}
	return Sum;
	}

static SCORE SPGapFreqs(const FCOUNT Freqs[])
	{
#if TRACE
	Log("Freqs=");
	for (unsigned i = 0; i < 4; ++i)
		if (Freqs[i] != 0)
			Log(" %s=%.3g", GapTypeToStr(i), Freqs[i]);
	Log("\n");
#endif

	SCORE TotalOffDiag = 0;
	SCORE TotalDiag = 0;
	for (unsigned i = 0; i < 4; ++i)
		{
		const FCOUNT fi = Freqs[i];
		if (0 == fi)
			continue;
		const float *Row = GapScoreMatrix[i];
		SCORE diagt = fi*fi*Row[i];
		TotalDiag += diagt;
#if	TRACE
		Log("SPFGaps %s %s + Mx=%.3g TotalDiag += %.3g\n",
		  GapTypeToStr(i),
		  GapTypeToStr(i),
		  Row[i],
		  diagt);
#endif
		SCORE Sum = 0;
		for (unsigned j = 0; j < i; ++j)
			{
			SCORE t = Freqs[j]*Row[j];
#if	TRACE
			if (Freqs[j] != 0)
				Log("SPFGaps %s %s + Mx=%.3g Sum += %.3g\n",
				  GapTypeToStr(i),
				  GapTypeToStr(j),
				  Row[j],
				  fi*t);
#endif
			Sum += t;
			}
		TotalOffDiag += fi*Sum;
		}
#if TRACE
	Log("SPFGap TotalOffDiag=%.3g + TotalDiag=%.3g = %.3g\n",
	  TotalOffDiag, TotalDiag, TotalOffDiag + TotalDiag);
#endif
	return TotalOffDiag*2 + TotalDiag;
	}

static SCORE SPFreqs(const FCOUNT Freqs[])
	{
#if TRACE
	Log("Freqs=");
	for (unsigned i = 0; i < 20; ++i)
		if (Freqs[i] != 0)
			Log(" %c=%.3g", LetterToCharAmino(i), Freqs[i]);
	Log("\n");
#endif

	SCORE TotalOffDiag = 0;
	SCORE TotalDiag = 0;
	for (unsigned i = 0; i < 20; ++i)
		{
		const FCOUNT fi = Freqs[i];
		if (0 == fi)
			continue;
		const float *Row = (*g_ptrScoreMatrix)[i];
		SCORE diagt = fi*fi*Row[i];
		TotalDiag += diagt;
#if	TRACE
		Log("SPF %c %c + Mx=%.3g TotalDiag += %.3g\n",
		  LetterToCharAmino(i),
		  LetterToCharAmino(i),
		  Row[i],
		  diagt);
#endif
		SCORE Sum = 0;
		for (unsigned j = 0; j < i; ++j)
			{
			SCORE t = Freqs[j]*Row[j];
#if	TRACE
			if (Freqs[j] != 0)
				Log("SPF %c %c + Mx=%.3g Sum += %.3g\n",
				  LetterToCharAmino(i),
				  LetterToCharAmino(j),
				  Row[j],
				  fi*t);
#endif
			Sum += t;
			}
		TotalOffDiag += fi*Sum;
		}
#if TRACE
	Log("SPF TotalOffDiag=%.3g + TotalDiag=%.3g = %.3g\n",
	  TotalOffDiag, TotalDiag, TotalOffDiag + TotalDiag);
#endif
	return TotalOffDiag*2 + TotalDiag;
	}

static SCORE ObjScoreSPCol(const MSA &msa, unsigned uColIndex)
	{
	FCOUNT Freqs[20];
	FCOUNT GapFreqs[4];

	memset(Freqs, 0, sizeof(Freqs));
	memset(GapFreqs, 0, sizeof(GapFreqs));

	const unsigned uSeqCount = msa.GetSeqCount();
#if	TRACE
	Log("Weights=");
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		Log(" %u=%.3g", uSeqIndex, msa.GetSeqWeight(uSeqIndex));
	Log("\n");
#endif
	SCORE SelfOverCount = 0;
	SCORE GapSelfOverCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		WEIGHT w = msa.GetSeqWeight(uSeqIndex);

		bool bGapThisCol = msa.IsGap(uSeqIndex, uColIndex);
		bool bGapPrevCol = (uColIndex == 0 ? false : msa.IsGap(uSeqIndex, uColIndex - 1));
		int GapType = bGapThisCol + 2*bGapPrevCol;
		assert(GapType >= 0 && GapType < 4);
		GapFreqs[GapType] += w;
		SCORE gapt = w*w*GapScoreMatrix[GapType][GapType];
		GapSelfOverCount += gapt;

		if (bGapThisCol)
			continue;
		unsigned uLetter = msa.GetLetterEx(uSeqIndex, uColIndex);
		if (uLetter >= 20)
			continue;
		Freqs[uLetter] += w;
		SCORE t = w*w*(*g_ptrScoreMatrix)[uLetter][uLetter];
#if	TRACE
		Log("FastCol compute freqs & SelfOverCount %c w=%.3g M=%.3g SelfOverCount += %.3g\n",
		  LetterToCharAmino(uLetter), w, (*g_ptrScoreMatrix)[uLetter][uLetter], t);
#endif
		SelfOverCount += t;
		}
	SCORE SPF = SPFreqs(Freqs);
	SCORE Col = SPF - SelfOverCount;

	SCORE SPFGaps = SPGapFreqs(GapFreqs);
	SCORE ColGaps = SPFGaps - GapSelfOverCount;
#if	TRACE
	Log("SPF=%.3g - SelfOverCount=%.3g = %.3g\n", SPF, SelfOverCount, Col);
	Log("SPFGaps=%.3g - GapsSelfOverCount=%.3g = %.3g\n", SPFGaps, GapSelfOverCount, ColGaps);
#endif
	return Col + ColGaps;
	}

SCORE ObjScoreSPDimer(const MSA &msa)
	{
	static bool bGapScoreMatrixInit = false;
	if (!bGapScoreMatrixInit)
		InitGapScoreMatrix();

	SCORE Total = 0;
	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		SCORE Col = ObjScoreSPCol(msa, uColIndex);
#if	TRACE
		{
		SCORE ColCheck = SPColBrute(msa, uColIndex);
		Log("FastCol=%.3g CheckCol=%.3g\n", Col, ColCheck);
		}
#endif
		Total += Col;
		}
#if TRACE
	Log("Total/2 = %.3g (final result from fast)\n", Total/2);
#endif
	return Total/2;
	}
