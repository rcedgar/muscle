#include "muscle.h"
#include "distfunc.h"
#include "msa.h"
#include "seqvect.h"
#include "pwpath.h"

void DistPWKimura(const SeqVect &v, DistFunc &DF)
	{
	SEQWEIGHT SeqWeightSave = GetSeqWeightMethod();
	SetSeqWeightMethod(SEQWEIGHT_Henikoff);

	const unsigned uSeqCount = v.Length();
	DF.SetCount(uSeqCount);

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	SetProgressDesc("PWKimura distance");
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const Seq &s1 = v.GetSeq(uSeqIndex1);
		MSA msa1;
		msa1.FromSeq(s1);
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
			{
			if (0 == uCount%20)
				Progress(uCount, uPairCount);
			++uCount;
			const Seq &s2 = v.GetSeq(uSeqIndex2);
			MSA msa2;
			msa2.FromSeq(s2);
		
			PWPath Path;
			MSA msaOut;
			AlignTwoMSAs(msa1, msa2, msaOut, Path, false, false);

			double dPctId = msaOut.GetPctIdentityPair(0, 1);
			float f = (float) KimuraDist(dPctId);

			DF.SetDist(uSeqIndex1, uSeqIndex2, f);
			}
		}
	ProgressStepsDone();

	SetSeqWeightMethod(SeqWeightSave);
	}
