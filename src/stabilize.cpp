#include "muscle.h"
#include "msa.h"

void Stabilize(const MSA &msa, MSA &msaStable)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();

	msaStable.SetSize(uSeqCount, uColCount);
	for (unsigned uId = 0; uId < uSeqCount; ++uId)
		{
		const unsigned uSeqIndex = msa.GetSeqIndex(uId);
		msaStable.SetSeqName(uId, msa.GetSeqName(uSeqIndex));
		msaStable.SetSeqId(uSeqIndex, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msa.GetChar(uSeqIndex, uColIndex);
			msaStable.SetChar(uId, uColIndex, c);
			}
		}
	}
