#include "muscle.h"
#include "upgma5.h"
#include "textfile.h"
#include "tree.h"
#include "omplock.h"

double GetProtDist(const char *Q, const char *T, uint ColCount);

static void MakeDistMx(const MultiSequence &Aln,
  vector<vector<float> > &DistMx, vector<string> &Labels)
	{
	DistMx.clear();
	Labels.clear();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();

	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount, FLT_MAX);
		DistMx[i][i] = 0;

		const string &Label = Aln.GetLabel(i);
		Labels.push_back(Label);
		}

	const uint ThreadCount = GetRequestedThreadCount();
	const uint PairCount = (SeqCount*(SeqCount-1))/2;
	uint SeqIndexi = UINT_MAX;
	uint SeqIndexj = UINT_MAX;
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		uint ThreadLocali = UINT_MAX;
		uint ThreadLocalj = UINT_MAX;

		LOCK();
		if (SeqIndexi == UINT_MAX)
			{
			SeqIndexi = 1;
			SeqIndexj = 0;
			}
		else
			{
			++SeqIndexj;
			if (SeqIndexj == SeqIndexi)
				{
				++SeqIndexi;
				SeqIndexj = 0;
				}
			}
		ThreadLocali = SeqIndexi;
		ThreadLocalj = SeqIndexj;
		ProgressStep(PairCounter++, PairCount, "Protdists");
		UNLOCK();

		uint Li;
		const char *Seqi = (const char *) Aln.GetByteSeq(ThreadLocali, Li);

		uint Lj;
		const char *Seqj = (const char *) Aln.GetByteSeq(ThreadLocalj, Lj);

		float dij = (float) GetProtDist(Seqi, Seqj, ColCount);

		DistMx[ThreadLocali][ThreadLocalj] = dij;
		DistMx[ThreadLocalj][ThreadLocali] = dij;
		}
	}

void cmd_upgma5_msa()
	{
	const string &InputFileName = opt(upgma5_msa);
	const string &OutputFileName = opt(output);

	MultiSequence Aln;
	Aln.FromFASTA(InputFileName);
	bool IsNucleo = Aln.GuessIsNucleo();
	SetAlphab(IsNucleo);

	UPGMA5 U;
	vector<vector<float> > &DistMx = U.m_DistMx;
	vector<string> &Labels = U.m_Labels;

	MakeDistMx(Aln, DistMx, Labels);

	LINKAGE Linkage = LINKAGE_Avg;
	string sLink = "avg";
	if (optset_linkage)
		{
		sLink = opt(linkage);
		if (sLink == "avg")
			Linkage = LINKAGE_Avg;
		else if (sLink == "min")
			Linkage = LINKAGE_Min;
		else if (sLink == "max")
			Linkage = LINKAGE_Max;
		else if (sLink == "biased")
			Linkage = LINKAGE_Biased;
		else
			Die("Invalid -linkage %s", sLink.c_str());
		}
	ProgressLog("UPGMA5(%s)\n", sLink.c_str());

	Tree t;
	U.Run(Linkage, t);
	TextFile TF(OutputFileName.c_str(), true);
	t.ToFile(TF);
	TF.Close();

	ProgressLog("All done.\n");
	}
