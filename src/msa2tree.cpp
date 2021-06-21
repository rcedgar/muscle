#include "muscle.h"
#include "msa.h"
#include "distfunc.h"
#include "distcalc.h"
#include "msa.h"
#include "tree.h"
#include "upgma5.h"

void cmd_msa2tree()
	{
	const string &MSAFileName = opt(msa2tree);
	const string &OutputFileName = opt(output);

	DISTANCE Dist = DISTANCE_PctIdKimura;
	if (optset_distance)
		{
		const string &sd = opt(distance);
		if (sd == "pctidkimura")
			Dist = DISTANCE_PctIdKimura;
		else if (sd == "pctidlog")
			Dist = DISTANCE_PctIdLog;
		else if (sd == "mydist")
			Dist = DISTANCE_MyDist;
		else if (sd == "kmer6_6")
			Dist = DISTANCE_Kmer6_6;
		else if (sd == "kmer20_3")
			Dist = DISTANCE_Kmer20_3;
		else
			Die("Invalid distance '%s'", sd.c_str());
		}

	LINKAGE Linkage = LINKAGE_Biased;
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

	const string DistStr = string(DISTANCEToStr(Dist));
	const string LinkStr = string(LINKAGEToStr(Linkage));

	ProgressLog("distance = %s\n", DistStr.c_str());
	ProgressLog("linkage = %s\n", LinkStr.c_str());

	MSA msa;
	Progress("Reading %s...", MSAFileName.c_str());
	msa.FromFASTAFile(MSAFileName);
	const uint SeqCount = msa.GetSeqCount();
	Progress(" %u seqs\n", SeqCount);

	ALPHA alpha = msa.GuessAlpha();
	SetAlpha(alpha);

	vector<string> Labels;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Name = msa.GetSeqName(SeqIndex);
		string Label = string(Name);
		Labels.push_back(Label);
		}

	DistCalcMSA DC;
	DC.Init(msa, Dist);

	vector<vector<float> > DistMx(SeqCount);
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		DistMx[SeqIndex1].resize(SeqCount, 0);

	const uint PairCount = (SeqCount*(SeqCount-1))/2;
	uint PairIndex = 0;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		for (uint SeqIndex2 = 0; SeqIndex2 < SeqIndex1; ++SeqIndex2)
			{
			ProgressStep(PairIndex++, PairCount, "Dist mx %s", DistStr.c_str());
			float d = DC.CalcDist(SeqIndex1, SeqIndex2);
			DistMx[SeqIndex1][SeqIndex2] = d;
			DistMx[SeqIndex2][SeqIndex1] = d;
			}
		}

	UPGMA5 U;
	U.Init(Labels, DistMx);
	U.ScaleDistMx();

	Tree tree;
	U.Run(Linkage, tree);

	tree.ToFile(OutputFileName);
	}
