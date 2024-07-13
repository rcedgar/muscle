#include "muscle.h"
#include "upgma5.h"
#include "pprog.h"
#include "super7.h"
#include "mpcflat_mega.h"
#include "mega.h"

void CalcGuideTree_SW_BLOSUM62(const MultiSequence &Input, Tree &T);

void cmd_super7_mega()
	{
	uint ShrubSize = 32;
	if (optset_shrub_size)
		ShrubSize = opt(shrub_size);
	if (ShrubSize < 3)
		Die("-shrub_size must be >= 3");

	Mega MM;
	MM.FromFile(g_Arg1);

	MultiSequence InputSeqs;
	InputSeqs.FromStrings(MM.m_Labels, MM.m_Seqs);
	SetGlobalInputMS(InputSeqs);

	bool Nucleo = InputSeqs.GuessIsNucleo();
	CheckMegaOpts(Nucleo);

	Tree GuideTree;
	if (optset_guidetreein)
		GuideTree.FromFile(opt(guidetreein));
	else if (optset_distmxin)
		{
		UPGMA5 U;
		U.ReadDistMx2(opt(distmxin));
		U.Run(LINKAGE_Avg, GuideTree);
		}
	else
		CalcGuideTree_SW_BLOSUM62(InputSeqs, GuideTree);

	SetAlpha(ALPHA_Amino);
	InitProbcons();

	Super7 S7;
	MPCFlat_mega *MPCm = new MPCFlat_mega;
	MPCm->m_MM = &MM;
	S7.m_MPC = MPCm;
	S7.Run(InputSeqs, GuideTree, ShrubSize);
	S7.m_FinalMSA.ToFasta(opt(output));

	ProgressLog("Done.\n");
	}
