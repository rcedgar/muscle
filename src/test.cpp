#include "muscle.h"
#include "m3alnparams.h"
#include "simplecluster.h"
#if 0
// Example from
// https://en.wikipedia.org/wiki/UPGMA
void cmd_test()
	{
	vector<string> Labels;
	vector<vector<float> > DistMx(5);

	uint i = 0;
#define row(Label, d0, d1, d2, d3, d4)	\
	Labels.push_back(#Label);			\
	DistMx[i].resize(5, FLT_MAX);		\
	DistMx[i][0] = d0;				\
	DistMx[i][1] = d1;				\
	DistMx[i][2] = d2;				\
	DistMx[i][3] = d3;				\
	DistMx[i][4] = d4;				\
	++i;								

	row(a,  0, 17, 21, 31, 23)
	row(b, 17,  0, 30, 34, 21)
	row(c, 21, 30,  0, 28, 39)
	row(d, 31, 34, 28,  0, 43)
	row(e, 23, 21, 39, 43,  0)
#undef row

	SimpleCluster SC;
	SC.Run(DistMx, Labels, "avg", false);

	Tree T;
	SC.GetTree(T);

	T.LogMe();
	T.ToFile("test.newick");
	}
#endif // 0
