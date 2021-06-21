#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "upgma5.h"
#include "distcalc.h"

void TreeFromMSA5(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root)
	{
	LINKAGE Linkage = LINKAGE_Undefined;
	switch (Cluster)
		{
	case CLUSTER_UPGMA:
		Linkage = LINKAGE_Avg;
		break;
	case CLUSTER_UPGMAMin:
		Linkage = LINKAGE_Min;
		break;
	case CLUSTER_UPGMAMax:
		Linkage = LINKAGE_Max;
		break;
	case CLUSTER_UPGMB:
		Linkage = LINKAGE_Biased;
		break;
	default:
		Die("TreeFromMSA_UPGMA, CLUSTER_%u not supported", Cluster);
		}

	asserta(Root == ROOT_FromClustering);
	const uint SeqCount = msa.GetSeqCount();
	vector<string> Labels;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Name = msa.GetSeqName(SeqIndex);
		string Label = string(Name);
		Labels.push_back(Label);
		}

	DistCalcMSA DC;
	DC.Init(msa, Distance);

	vector<vector<float> > DistMx(SeqCount);
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		DistMx[SeqIndex1].resize(SeqCount, 0);

	const string DistStr = string(DISTANCEToStr(Distance));
	const string LinkStr = string(LINKAGEToStr(Linkage));

	const uint PairCount = (SeqCount*(SeqCount-1))/2;
	uint PairIndex = 0;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		for (uint SeqIndex2 = 0; SeqIndex2 < SeqIndex1; ++SeqIndex2)
			{
			ProgressStep(PairIndex++, PairCount, "Dist mx %s/%s",
			  DistStr.c_str(), LinkStr.c_str());
			float d = DC.CalcDist(SeqIndex1, SeqIndex2);
			DistMx[SeqIndex1][SeqIndex2] = d;
			DistMx[SeqIndex2][SeqIndex1] = d;
			}
		}

	UPGMA5 U;
	U.Init(Labels, DistMx);
	U.ScaleDistMx();

	U.Run(Linkage, tree);
	}
