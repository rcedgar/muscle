#include "muscle.h"
#include "msa.h"

void OutWeights(const char *FileName, const MSA &msa)
	{
	FILE *f = fopen(FileName, "w");
	if (0 == f)
		Quit("Cannot open '%s'", FileName);
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const char *Id = msa.GetSeqName(uSeqIndex);
		const WEIGHT w = msa.GetSeqWeight(uSeqIndex);
		fprintf(f, "%s\t%.3g\n", Id, w);
		}
	fclose(f);
	}
