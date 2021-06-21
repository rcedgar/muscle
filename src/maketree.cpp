#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "tree.h"

void DoMakeTree()
	{
	if (g_pstrInFileName == 0 || g_pstrOutFileName == 0)
		Quit("-maketree requires -in <msa> and -out <treefile>");

	SetStartTime();

	SetSeqWeightMethod(g_SeqWeight1);

	TextFile MSAFile(g_pstrInFileName);

	MSA msa;
	msa.FromFile(MSAFile);

	unsigned uSeqCount = msa.GetSeqCount();
	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	SetMuscleInputMSA(msa);

	Progress("%u sequences", uSeqCount);

	Tree tree;
	TreeFromMSA(msa, tree, g_Cluster2, g_Distance2, g_Root2);

	TextFile TreeFile(g_pstrOutFileName, true);
	tree.ToFile(TreeFile);

	Progress("Tree created");
	}
