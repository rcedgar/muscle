#include "muscle.h"
#include "uclustpd.h"
#include "upgma5.h"

void cmd_build_guide_tree()
	{
#if 0
	const string &InputFileName = opt(build_guide_tree);
	const string &TreeFileName = opt(output);
	asserta(optset_maxpd);
	double MaxPD = opt(maxpd);

	MultiSequence Input;
	Input.FromFASTA(InputFileName, true);
	const uint SeqCount = Input.GetSeqCount();

	bool IsNucleo = Input.GuessIsNucleo();
	SetSubstMx(IsNucleo);

	vector<uint> AllSeqIndexes;
	for (uint i = 0; i < SeqCount; ++i)
		AllSeqIndexes.push_back(i);

	UClustPD UD;
	UD.Run(Input, AllSeqIndexes, MaxPD);

	//FILE *fOut = CreateStdioFile(opt(tsvout));
	//UD.ToTsv(fOut);
	//CloseStdioFile(fOut);

	uint MaxClusterSize = 512;

	const uint ClusterCount = UD.GetClusterCount();
	vector<UClustPD *> UCs(ClusterCount);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount;
	  ++ClusterIndex)
		{
		uint ClusterSize = UD.GetClusterSize(ClusterIndex);
		if (ClusterSize <= MaxClusterSize)
			{
			UClustPD *SubUC = new UClustPD;
			vector<uint> ClusterSeqIndexes;
			UD.GetClusterSeqIndexes(ClusterIndex, ClusterSeqIndexes);
			SubUC->Run(Input, ClusterSeqIndexes, MaxPD/2);
			}
		else
			{
			}
		}
#endif // 0
	}
