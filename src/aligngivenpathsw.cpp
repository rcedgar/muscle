#include "muscle.h"
#include "msa.h"
#include "pwpath.h"
#include "profile.h"

#define	TRACE	0

static void AppendDelete(const MSA &msaA, unsigned &uColIndexA,
  unsigned uSeqCountA, unsigned uSeqCountB, MSA &msaCombined,
  unsigned &uColIndexCombined)
	{
#if	TRACE
	Log("AppendDelete ColIxA=%u ColIxCmb=%u\n",
	  uColIndexA, uColIndexCombined);
#endif
	for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
		{
		char c = msaA.GetChar(uSeqIndexA, uColIndexA);
		msaCombined.SetChar(uSeqIndexA, uColIndexCombined, c);
		}

	for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
		msaCombined.SetChar(uSeqCountA + uSeqIndexB, uColIndexCombined, '-');

	++uColIndexCombined;
	++uColIndexA;
	}

static void AppendInsert(const MSA &msaB, unsigned &uColIndexB,
  unsigned uSeqCountA, unsigned uSeqCountB, MSA &msaCombined,
  unsigned &uColIndexCombined)
	{
#if	TRACE
	Log("AppendInsert ColIxB=%u ColIxCmb=%u\n",
	  uColIndexB, uColIndexCombined);
#endif
	for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
		msaCombined.SetChar(uSeqIndexA, uColIndexCombined, '-');

	for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
		{
		char c = msaB.GetChar(uSeqIndexB, uColIndexB);
		msaCombined.SetChar(uSeqCountA + uSeqIndexB, uColIndexCombined, c);
		}

	++uColIndexCombined;
	++uColIndexB;
	}

static void AppendUnalignedTerminals(const MSA &msaA, unsigned &uColIndexA, unsigned uColCountA,
  const MSA &msaB, unsigned &uColIndexB, unsigned uColCountB, unsigned uSeqCountA,
  unsigned uSeqCountB, MSA &msaCombined, unsigned &uColIndexCombined)
	{
#if	TRACE
	Log("AppendUnalignedTerminals ColIxA=%u ColIxB=%u ColIxCmb=%u\n",
	  uColIndexA, uColIndexB, uColIndexCombined);
#endif
	const unsigned uLengthA = msaA.GetColCount();
	const unsigned uLengthB = msaB.GetColCount();

	unsigned uNewColCount = uColCountA;
	if (uColCountB > uNewColCount)
		uNewColCount = uColCountB;

	for (unsigned n = 0; n < uColCountA; ++n)
		{
		for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
			{
			char c = msaA.GetChar(uSeqIndexA, uColIndexA + n);
			c = UnalignChar(c);
			msaCombined.SetChar(uSeqIndexA, uColIndexCombined + n, c);
			}
		}
	for (unsigned n = uColCountA; n < uNewColCount; ++n)
		{
		for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
			msaCombined.SetChar(uSeqIndexA, uColIndexCombined + n, '.');
		}

	for (unsigned n = 0; n < uColCountB; ++n)
		{
		for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
			{
			char c = msaB.GetChar(uSeqIndexB, uColIndexB + n);
			c = UnalignChar(c);
			msaCombined.SetChar(uSeqCountA + uSeqIndexB, uColIndexCombined + n, c);
			}
		}
	for (unsigned n = uColCountB; n < uNewColCount; ++n)
		{
		for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
			msaCombined.SetChar(uSeqCountA + uSeqIndexB, uColIndexCombined + n, '.');
		}

	uColIndexCombined += uNewColCount;
	uColIndexA += uColCountA;
	uColIndexB += uColCountB;
	}

static void AppendMatch(const MSA &msaA, unsigned &uColIndexA, const MSA &msaB,
  unsigned &uColIndexB, unsigned uSeqCountA, unsigned uSeqCountB,
  MSA &msaCombined, unsigned &uColIndexCombined)
	{
#if	TRACE
	Log("AppendMatch ColIxA=%u ColIxB=%u ColIxCmb=%u\n",
	  uColIndexA, uColIndexB, uColIndexCombined);
#endif

	for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
		{
		char c = msaA.GetChar(uSeqIndexA, uColIndexA);
		msaCombined.SetChar(uSeqIndexA, uColIndexCombined, c);
		}

	for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
		{
		char c = msaB.GetChar(uSeqIndexB, uColIndexB);
		msaCombined.SetChar(uSeqCountA + uSeqIndexB, uColIndexCombined, c);
		}

	++uColIndexA;
	++uColIndexB;
	++uColIndexCombined;
	}

void AlignTwoMSAsGivenPathSW(const PWPath &Path, const MSA &msaA, const MSA &msaB,
  MSA &msaCombined)
	{
	msaCombined.Clear();

#if	TRACE
	Log("AlignTwoMSAsGivenPathSW\n");
	Log("Template A:\n");
	msaA.LogMe();
	Log("Template B:\n");
	msaB.LogMe();
#endif

	const unsigned uColCountA = msaA.GetColCount();
	const unsigned uColCountB = msaB.GetColCount();

	const unsigned uSeqCountA = msaA.GetSeqCount();
	const unsigned uSeqCountB = msaB.GetSeqCount();

	msaCombined.SetSeqCount(uSeqCountA + uSeqCountB);

// Copy sequence names into combined MSA
	for (unsigned uSeqIndexA = 0; uSeqIndexA < uSeqCountA; ++uSeqIndexA)
		{
		msaCombined.SetSeqName(uSeqIndexA, msaA.GetSeqName(uSeqIndexA));
		msaCombined.SetSeqId(uSeqIndexA, msaA.GetSeqId(uSeqIndexA));
		}

	for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
		{
		msaCombined.SetSeqName(uSeqCountA + uSeqIndexB, msaB.GetSeqName(uSeqIndexB));
		msaCombined.SetSeqId(uSeqCountA + uSeqIndexB, msaB.GetSeqId(uSeqIndexB));
		}

	unsigned uColIndexA = 0;
	unsigned uColIndexB = 0;
	unsigned uColIndexCombined = 0;
	const unsigned uEdgeCount = Path.GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
#if	TRACE
		Log("\nEdge %u %c%u.%u\n",
		  uEdgeIndex,
		  Edge.cType,
		  Edge.uPrefixLengthA,
		  Edge.uPrefixLengthB);
#endif
		const char cType = Edge.cType;
		const unsigned uPrefixLengthA = Edge.uPrefixLengthA;
		unsigned uColCountA = 0;
		if (uPrefixLengthA > 0)
			{
			const unsigned uNodeIndexA = uPrefixLengthA - 1;
			const unsigned uTplColIndexA = uNodeIndexA;
			if (uTplColIndexA > uColIndexA)
				uColCountA = uTplColIndexA - uColIndexA;
			}

		const unsigned uPrefixLengthB = Edge.uPrefixLengthB;
		unsigned uColCountB = 0;
		if (uPrefixLengthB > 0)
			{
			const unsigned uNodeIndexB = uPrefixLengthB - 1;
			const unsigned uTplColIndexB = uNodeIndexB;
			if (uTplColIndexB > uColIndexB)
				uColCountB = uTplColIndexB - uColIndexB;
			}

		AppendUnalignedTerminals(msaA, uColIndexA, uColCountA, msaB, uColIndexB, uColCountB,
		  uSeqCountA, uSeqCountB, msaCombined, uColIndexCombined);

		switch (cType)
			{
		case 'M':
			{
			assert(uPrefixLengthA > 0);
			assert(uPrefixLengthB > 0);
			const unsigned uColA = uPrefixLengthA - 1;
			const unsigned uColB = uPrefixLengthB - 1;
			assert(uColIndexA == uColA);
			assert(uColIndexB == uColB);
			AppendMatch(msaA, uColIndexA, msaB, uColIndexB, uSeqCountA, uSeqCountB,
			  msaCombined, uColIndexCombined);
			break;
			}
		case 'D':
			{
			assert(uPrefixLengthA > 0);
			const unsigned uColA = uPrefixLengthA - 1;
			assert(uColIndexA == uColA);
			AppendDelete(msaA, uColIndexA, uSeqCountA, uSeqCountB, msaCombined, uColIndexCombined);
			break;
			}
		case 'I':
			{
			assert(uPrefixLengthB > 0);
			const unsigned uColB = uPrefixLengthB - 1;
			assert(uColIndexB == uColB);
			AppendInsert(msaB, uColIndexB, uSeqCountA, uSeqCountB, msaCombined, uColIndexCombined);
			break;
			}
		default:
			assert(false);
			}
		}
	unsigned uInsertColCountA = uColCountA - uColIndexA;
	unsigned uInsertColCountB = uColCountB - uColIndexB;

	AppendUnalignedTerminals(msaA, uColIndexA, uInsertColCountA, msaB, uColIndexB,
	  uInsertColCountB, uSeqCountA, uSeqCountB, msaCombined, uColIndexCombined);
	}
