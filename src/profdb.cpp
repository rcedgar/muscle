#include "muscle.h"
#include "textfile.h"
#include "seqvect.h"
#include "distfunc.h"
#include "msa.h"
#include "tree.h"
#include "clust.h"
#include "profile.h"
#include "clustsetmsa.h"

void ProfDB()
	{
	SetOutputFileName(g_pstrOutFileName);
	SetInputFileName(g_pstrFileName2);
	SetStartTime();

	TextFile file1(g_pstrFileName1);
	TextFile file2(g_pstrFileName2);

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	TextFile fileIn(g_pstrFileName1);
	MSA msa1;
	msa1.FromFile(fileIn);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	if (0 == uSeqCount1)
		Quit("No sequences in input alignment");

	SeqVect v;
	v.FromFASTAFile(file2);
	const unsigned uSeqCount2 = v.Length();
	if (0 == uSeqCount2)
		Quit("No sequences in input alignment");

	MSA::SetIdCount(uSeqCount1 + uSeqCount2);
	SetProgressDesc("Align sequence database to profile");
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount2; ++uSeqIndex)
		{
		Progress(uSeqIndex, uSeqCount2);
		Seq &s = *(v[uSeqIndex]);
		s.SetId(0);
		MSA msaTmp;
		msaTmp.FromSeq(s);
		MSA msaOut;
		ProfileProfile(msa1, msaTmp, msaOut);
		msa1.Copy(msaOut);
		}
	ProgressStepsDone();

	TextFile fileOut(g_pstrOutFileName, true);
	msa1.ToFile(fileOut);
	}
