#include "myutils.h"
#include "muscle.h"
#include "seqvect.h"
#include "probcons.h"
#include "upgma5.h"

void DistMxFromSeqVect_EA(const SeqVect &SV, vector<vector<float> > &DistMx,
  vector<string> &Labels);

void TreeFromSeqVect_EA(const SeqVect &SV, Tree &Tree)
	{
	vector<string> Labels;
	vector<vector<float> > DistMx;
	DistMxFromSeqVect_EA(SV, DistMx, Labels);

	UPGMA5 U;
	U.Init(Labels, DistMx);
	U.Run(LINKAGE_Biased, Tree);
	}
