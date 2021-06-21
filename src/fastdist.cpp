#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"

void DistPWScoreDist(const SeqVect &v, DistFunc &DF);

void DistUnaligned(const SeqVect &v, DISTANCE DistMethod, DistFunc &DF)
	{
	const unsigned uSeqCount = v.Length();

	switch (DistMethod)
		{
	case DISTANCE_Kmer6_6:
		DistKmer6_6(v, DF);
		break;

	case DISTANCE_Kmer20_3:
		DistKmer20_3(v, DF);
		break;

	case DISTANCE_Kmer20_4:
		FastDistKmer(v, DF);
		break;

	case DISTANCE_Kbit20_3:
		DistKbit20_3(v, DF);
		break;

	case DISTANCE_Kmer4_6:
		DistKmer4_6(v, DF);
		break;

	case DISTANCE_PWKimura:
		DistPWKimura(v, DF);
		break;

	case DISTANCE_PWScoreDist:
		DistPWScoreDist(v, DF);
		break;

	default:
		Quit("DistUnaligned, unsupported distance method %d", DistMethod);
		}

//	const char **SeqNames = (const char **) malloc(uSeqCount*sizeof(char *));
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq &s = *(v[uSeqIndex]);

		const char *ptrName = s.GetName();
		unsigned uId = s.GetId();

		DF.SetName(uSeqIndex, ptrName);
		DF.SetId(uSeqIndex, uId);
		}
	}
