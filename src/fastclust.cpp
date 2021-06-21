#include "muscle.h"
#include "seqvect.h"
#include "distfunc.h"
#include "clust.h"
#include "clustsetdf.h"
#include "tree.h"
#include "clust.h"
#include "distcalc.h"
#include <math.h>

static void TreeFromSeqVect_UPGMA(const DistFunc &DF, CLUSTER Cluster, Tree &tree)
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
		Quit("TreeFromSeqVect_UPGMA, CLUSTER_%u not supported", Cluster);
		}
	
	DistCalcDF DC;
	DC.Init(DF);
	UPGMA2(DC, tree, Linkage);
	}

void TreeFromSeqVect(const SeqVect &v, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root)
	{
	DistFunc DF;
	DistUnaligned(v, Distance, DF);
	DF.m_DistType = DISTANCEToStr(Distance);

	TreeFromSeqVect_UPGMA(DF, Cluster, tree);
	FixRoot(tree, Root);
	}
