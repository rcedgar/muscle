#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "quarts.h"
#include "alpha3.h"

static double CalcScore(const char *Seqi, const char *Seqj, uint L)
	{
	extern int BLOSUM62[20][20];

	double Score = 0;
	bool InGap = false;
	const double GAP_OPEN = -11;
	const double GAP_EXT = -1;
	for (uint i = 0; i < L; ++i)
		{
		byte ci = (byte) Seqi[i];
		byte cj = (byte) Seqj[i];
		bool gapi = isgap(ci);
		bool gapj = isgap(cj);
		if (gapi && gapj)
			continue;
		if (gapi || gapj)
			{
			if (InGap)
				Score += GAP_EXT;
			else
				{
				InGap = true;
				Score += GAP_OPEN;
				}
			continue;
			}
		uint Letteri = g_CharToLetterAmino[ci];
		uint Letterj = g_CharToLetterAmino[cj];
		if (Letteri >= 20 || Letterj >= 20)
			continue;
		InGap = false;
		int Sub = BLOSUM62[Letteri][Letterj];
		Score += Sub;
		}
	return Score;
	}

void cmd_blosumscore()
	{
	const string &MSAFileName = opt(blosumscore);
	MSA Aln;
	TextFile f(MSAFileName.c_str());
	Aln.FromFASTAFile(f);
	f.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();

	const uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	double Sum = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const char *Seqi = Aln.GetSeqCharPtr(i);
		for (uint j = i+1; j < SeqCount; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Calculating");
			const char *Seqj = Aln.GetSeqCharPtr(j);
			double Score = CalcScore(Seqi, Seqj, ColCount);
			Sum += Score;
			}
		}
	asserta(PairIndex == PairCount);

	ProgressPrefix(false);
	ProgressLog("\n");
	ProgressLog("Score = %.3g\n", Sum/PairCount);
	ProgressLog("\n");
	ProgressPrefix(true);
	}
