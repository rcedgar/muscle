#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "clust.h"
#include "clustsetmsa.h"
#include "distcalc.h"

static void TreeFromMSA_UPGMA(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance)
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
		Quit("TreeFromMSA_UPGMA, CLUSTER_%u not supported", Cluster);
		}
	
	DistCalcMSA DC;
	DC.Init(msa, Distance);
	UPGMA2(DC, tree, Linkage);
	}

void TreeFromMSA(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root)
	{
	TreeFromMSA_UPGMA(msa, tree, Cluster, Distance);
	FixRoot(tree, Root);
	}
