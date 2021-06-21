#include <limits.h>
#include <math.h>
#include "muscle.h"
#include "msa.h"
#include "distfunc.h"
#include "msa.h"
#include "seqvect.h"
#include "pwpath.h"

// ScoreDist
// E. Sonnhammer & V. Hollich, Scoredist: A simple and robust protein sequence
// distance estimator, BMC Bioinformatics 2005, 6:108.

extern int BLOSUM62[20][20];
extern double BLOSUM62_Expected;

static const double Dayhoff_CalibrationFactor = 1.3370;
static const double JTT_CalibrationFactor = 1.2873;
static const double MV_CalibrationFactor = 1.1775;
static const double LARGE_D = 3.0;

static double CalibrationFactor = JTT_CalibrationFactor;


// Similarity score
static double Sigma(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2,
  unsigned *ptrLength)
	{
	unsigned Length = 0;
	double Score = 0;
	const unsigned ColCount = msa.GetColCount();
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		unsigned Letter1 = msa.GetLetterEx(SeqIndex1, ColIndex);
		unsigned Letter2 = msa.GetLetterEx(SeqIndex2, ColIndex);
		if (Letter1 >= 20 || Letter2 >= 20)
			continue;
		++Length;
		Score += BLOSUM62[Letter1][Letter2];
		}

	*ptrLength = Length;
	return Score;
	}

// Normalized score
static double Sigma_N(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2)
	{
	unsigned Length = UINT_MAX;
	double Score = Sigma(msa, SeqIndex1, SeqIndex2, &Length);
	double RandomScore = Length*BLOSUM62_Expected;
	return Score - RandomScore;
	}

// Upper limit
static double Sigma_U(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2,
  unsigned *ptrLength)
	{
	double Score11 = Sigma(msa, SeqIndex1, SeqIndex1, ptrLength);
	double Score22 = Sigma(msa, SeqIndex2, SeqIndex2, ptrLength);
	return (Score11 + Score22)/2;
	}

// Normalized upper limit
static double Sigma_UN(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2)
	{
	unsigned Length = UINT_MAX;
	double Score = Sigma_U(msa, SeqIndex1, SeqIndex2, &Length);
	double RandomScore = Length*BLOSUM62_Expected;
	return Score - RandomScore;
	}

double GetScoreDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2)
	{
	if (g_Alpha != ALPHA_Amino)
		Quit("Scoredist is only for amino acid sequences");

	double s_N = Sigma_N(msa, SeqIndex1, SeqIndex2);
	double s_UN = Sigma_UN(msa, SeqIndex1, SeqIndex2);
	double d = 0.0;
	if (s_UN != 0)
		{
		double Ratio = s_N/s_UN;
		if (Ratio < 0.001)
			d = LARGE_D;
		else
			d =	-log(Ratio);
		}
	return d*CalibrationFactor;
	}

void DistPWScoreDist(const SeqVect &v, DistFunc &DF)
	{
	SEQWEIGHT SeqWeightSave = GetSeqWeightMethod();
	SetSeqWeightMethod(SEQWEIGHT_Henikoff);

	const unsigned uSeqCount = v.Length();
	DF.SetCount(uSeqCount);

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	SetProgressDesc("PW ScoreDist");
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const Seq &s1 = v.GetSeq(uSeqIndex1);
		MSA msa1;
		msa1.FromSeq(s1);
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
			{
			if (0 == uCount%20)
				Progress(uCount, uPairCount);
			++uCount;
			const Seq &s2 = v.GetSeq(uSeqIndex2);
			MSA msa2;
			msa2.FromSeq(s2);
		
			PWPath Path;
			MSA msaOut;
			AlignTwoMSAs(msa1, msa2, msaOut, Path, false, false);

			float d = (float) GetScoreDist(msaOut, 0, 1);
			DF.SetDist(uSeqIndex1, uSeqIndex2, d);
			}
		}
	ProgressStepsDone();

	SetSeqWeightMethod(SeqWeightSave);
	}
