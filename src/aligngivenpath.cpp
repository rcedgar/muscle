#include "muscle.h"
#include "msa.h"
#include "pwpath.h"
#include "profile.h"

#define	TRACE	0

static void LogPP(const ProfPos &PP)
	{
	Log("ResidueGroup   %u\n", PP.m_uResidueGroup);
	Log("AllGaps      %d\n", PP.m_bAllGaps);
	Log("Occ          %.3g\n", PP.m_fOcc);
	Log("LL=%.3g LG=%.3g GL=%.3g GG=%.3g\n", PP.m_LL, PP.m_LG, PP.m_GL, PP.m_GG);
	Log("Freqs        ");
	for (unsigned i = 0; i < 20; ++i)
		if (PP.m_fcCounts[i] > 0)
			Log("%c=%.3g ", LetterToChar(i), PP.m_fcCounts[i]);
	Log("\n");
	}

static void AssertProfPosEq(const ProfPos *PA, const ProfPos *PB, unsigned i)
	{
	const ProfPos &PPA = PA[i];
	const ProfPos &PPB = PB[i];
#define	eq(x)	if (PPA.m_##x != PPB.m_##x) { LogPP(PPA); LogPP(PPB); Quit("AssertProfPosEq." #x); }
#define be(x)	if (!BTEq(PPA.m_##x, PPB.m_##x)) { LogPP(PPA); LogPP(PPB); Quit("AssertProfPosEq." #x); }
	eq(bAllGaps)
	eq(uResidueGroup)

	be(LL)
	be(LG)
	be(GL)
	be(GG)
	be(fOcc)
	be(scoreGapOpen)
	be(scoreGapClose)

	for (unsigned j = 0; j < 20; ++j)
		{
#define	eqj(x)	if (PPA.m_##x != PPB.m_##x) Quit("AssertProfPosEq j=%u " #x, j);
#define bej(x)	if (!BTEq(PPA.m_##x, PPB.m_##x)) Quit("AssertProfPosEq j=%u " #x, j);
		bej(fcCounts[j]);
//		eqj(uSortOrder[j]) // may differ due to ties, don't check?
		bej(AAScores[j])
#undef eqj
#undef bej
		}
#undef eq
#undef be
	}

void AssertProfsEq(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB)
	{
	if (uLengthA != uLengthB)
		Quit("AssertProfsEq: lengths differ %u %u", uLengthA, uLengthB);
	for (unsigned i = 0; i < uLengthB; ++i)
		AssertProfPosEq(PA, PB, i);
	}

#if	DEBUG
static void ValidateProf(const ProfPos *Prof, unsigned uLength)
	{
	for (unsigned i = 0; i < uLength; ++i)
		{
		const ProfPos &PP = Prof[i];

		FCOUNT s1 = PP.m_LL + PP.m_LG + PP.m_GL + PP.m_GG;
		assert(BTEq(s1, 1.0));

		if (i > 0)
			{
			const ProfPos &PPPrev = Prof[i-1];
			FCOUNT s2 = PPPrev.m_LL + PPPrev.m_GL;
			FCOUNT s3 = PP.m_LL + PP.m_LG;
			assert(BTEq(s2, s3));
			}
		if (i < uLength - 1)
			{
			const ProfPos &PPNext = Prof[i+1];
			FCOUNT s4 = PP.m_LL + PP.m_GL;
			FCOUNT s5 = PPNext.m_LL + PPNext.m_LG;
			assert(BTEq(s4, s5));
			}
		}
	}
#else
#define ValidateProf(Prof, Length)	/* empty */
#endif

static void ScoresFromFreqsPos(ProfPos *Prof, unsigned uLength, unsigned uPos)
	{
	ProfPos &PP = Prof[uPos];
	SortCounts(PP.m_fcCounts, PP.m_uSortOrder);
	PP.m_uResidueGroup = ResidueGroupFromFCounts(PP.m_fcCounts);

// "Occupancy"
	PP.m_fOcc = PP.m_LL + PP.m_GL;

// Frequency of gap-opens in this position (i)
// Gap open 	= letter in i-1 and gap in i
//				= iff LG in i
	FCOUNT fcOpen = PP.m_LG;

// Frequency of gap-closes in this position
// Gap close	= gap in i and letter in i+1
//				= iff GL in i+1
	FCOUNT fcClose;
	if (uPos + 1 < uLength)
		fcClose = Prof[uPos + 1].m_GL;
	else
		fcClose = PP.m_GG + PP.m_LG;

	PP.m_scoreGapOpen = (SCORE) ((1.0 - fcOpen)*g_scoreGapOpen/2.0);
	PP.m_scoreGapClose = (SCORE) ((1.0 - fcClose)*g_scoreGapOpen/2.0);
#if	DOUBLE_AFFINE
	PP.m_scoreGapOpen2 = (SCORE) ((1.0 - fcOpen)*g_scoreGapOpen2/2.0);
	PP.m_scoreGapClose2 = (SCORE) ((1.0 - fcClose)*g_scoreGapOpen2/2.0);
#endif

	for (unsigned i = 0; i < g_AlphaSize; ++i)
		{
		SCORE scoreSum = 0;
		for (unsigned j = 0; j < g_AlphaSize; ++j)
			scoreSum += PP.m_fcCounts[j]*(*g_ptrScoreMatrix)[i][j];
		PP.m_AAScores[i] = scoreSum;
		}
	}

void ProfScoresFromFreqs(ProfPos *Prof, unsigned uLength)
	{
	for (unsigned i = 0; i < uLength; ++i)
		ScoresFromFreqsPos(Prof, uLength, i);
	}

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

static void AppendTplInserts(const MSA &msaA, unsigned &uColIndexA, unsigned uColCountA,
  const MSA &msaB, unsigned &uColIndexB, unsigned uColCountB, unsigned uSeqCountA,
  unsigned uSeqCountB, MSA &msaCombined, unsigned &uColIndexCombined)
	{
#if	TRACE
	Log("AppendTplInserts ColIxA=%u ColIxB=%u ColIxCmb=%u\n",
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

void AlignTwoMSAsGivenPath(const PWPath &Path, const MSA &msaA, const MSA &msaB,
  MSA &msaCombined)
	{
	msaCombined.Clear();

#if	TRACE
	Log("FastAlignProfiles\n");
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
		if (msaA.m_IdToSeqIndex != 0)
			msaCombined.SetSeqId(uSeqIndexA, msaA.GetSeqId(uSeqIndexA));
		}

	for (unsigned uSeqIndexB = 0; uSeqIndexB < uSeqCountB; ++uSeqIndexB)
		{
		msaCombined.SetSeqName(uSeqCountA + uSeqIndexB, msaB.GetSeqName(uSeqIndexB));
		if (msaB.m_IdToSeqIndex != 0)
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

// TODO: This code looks like a hangover from HMM estimation -- can we delete it?
		assert(uColCountA == 0);
		assert(uColCountB == 0);
		AppendTplInserts(msaA, uColIndexA, uColCountA, msaB, uColIndexB, uColCountB,
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

// TODO: This code looks like a hangover from HMM estimation -- can we delete it?
	assert(uInsertColCountA == 0);
	assert(uInsertColCountB == 0);
	AppendTplInserts(msaA, uColIndexA, uInsertColCountA, msaB, uColIndexB,
	  uInsertColCountB, uSeqCountA, uSeqCountB, msaCombined, uColIndexCombined);

	assert(msaCombined.GetColCount() == uEdgeCount);
	}

static const ProfPos PPStart =
	{
	false,		//m_bAllGaps;
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // m_uSortOrder[21];
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // m_fcCounts[20];
	1.0,	// m_LL;
	0.0,	// m_LG;
	0.0,	// m_GL;
	0.0,	// m_GG;
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // m_ALScores
	0,		// m_uResidueGroup;
	1.0,	// m_fOcc;
	0.0,	// m_fcStartOcc;
	0.0,	// m_fcEndOcc;
	0.0,	// m_scoreGapOpen;
	0.0,	// m_scoreGapClose;
	};

// MM
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetGapsMM(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wA*PPA.m_LL + wB*PPB.m_LL;
	PPO.m_LG = wA*PPA.m_LG + wB*PPB.m_LG;
	PPO.m_GL = wA*PPA.m_GL + wB*PPB.m_GL;
	PPO.m_GG = wA*PPA.m_GG + wB*PPB.m_GG;
	}

// MD
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  Bj		(-)
//  X		-	?L	LG
//  -		-	?G	GG
static void SetGapsMD(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wA*PPA.m_LL;
	PPO.m_LG = wA*PPA.m_LG + wB*(PPB.m_LL + PPB.m_GL);
	PPO.m_GL = wA*PPA.m_GL;
	PPO.m_GG = wA*PPA.m_GG + wB*(PPB.m_LG + PPB.m_GG);
	}

// DD
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  (-)		(-)
//  -		-	??	GG
static void SetGapsDD(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wA*PPA.m_LL;
	PPO.m_LG = wA*PPA.m_LG;
	PPO.m_GL = wA*PPA.m_GL;
	PPO.m_GG = wA*PPA.m_GG + wB;
	}

// MI
//  Ai		(-)		Out
//  X		-	?L	LG
//  -		-	?G	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetGapsMI(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wB*PPB.m_LL;
	PPO.m_LG = wB*PPB.m_LG + wA*(PPA.m_LL + PPA.m_GL);
	PPO.m_GL = wB*PPB.m_GL;
	PPO.m_GG = wB*PPB.m_GG + wA*(PPA.m_LG + PPA.m_GG);
	}

// DM
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  (-)		Bj		
//  -		X		?L	GL
//  -		-		?G	GG
static void SetGapsDM(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wA*PPA.m_LL;
	PPO.m_LG = wA*PPA.m_LG;
	PPO.m_GL = wA*PPA.m_GL + wB*(PPB.m_LL + PPB.m_GL);
	PPO.m_GG = wA*PPA.m_GG + wB*(PPB.m_LG + PPB.m_GG);
	}

// IM
//  (-)		Ai		Out		
//  -		X	?L	GL
//  -		-	?G	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetGapsIM(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wB*PPB.m_LL;
	PPO.m_LG = wB*PPB.m_LG;
	PPO.m_GL = wB*PPB.m_GL + wA*(PPA.m_LL + PPA.m_GL);
	PPO.m_GG = wB*PPB.m_GG + wA*(PPA.m_LG + PPA.m_GG);
	}

// ID
//  (-)		Ai		Out
//  -		X	?L	GL
//  -		-	?G	GG

//  Bj		(-)
//  X		-	?L	LG
//  -		-	?G	GG
static void SetGapsID(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = 0;
	PPO.m_LG = wB*PPB.m_GL + wB*PPB.m_LL;
	PPO.m_GL = wA*PPA.m_GL + wA*PPA.m_LL;
	PPO.m_GG = wA*(PPA.m_LG + PPA.m_GG) + wB*(PPB.m_LG + PPB.m_GG);
	}

// DI
//  Ai		(-)		Out
//  X		-	?L	LG
//  -		-	?G	GG

//  (-)		Bj
//  -		X	?L	GL
//  -		-	?G	GG
static void SetGapsDI(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = 0;
	PPO.m_LG = wA*PPA.m_GL + wA*PPA.m_LL;
	PPO.m_GL = wB*PPB.m_GL + wB*PPB.m_LL;
	PPO.m_GG = wA*(PPA.m_LG + PPA.m_GG) + wB*(PPB.m_LG + PPB.m_GG);
	}

// II
//  (-)		(-)		Out
//  -		-	??	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetGapsII(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	PPO.m_LL = wB*PPB.m_LL;
	PPO.m_LG = wB*PPB.m_LG;
	PPO.m_GL = wB*PPB.m_GL;
	PPO.m_GG = wB*PPB.m_GG + wA;
	}

static void SetFreqs(
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos *POut, unsigned uColIndexOut)
	{
	const ProfPos &PPA = uPrefixLengthA > 0 ? PA[uPrefixLengthA-1] : PPStart;
	const ProfPos &PPB = uPrefixLengthB > 0 ? PB[uPrefixLengthB-1] : PPStart;
	ProfPos &PPO = POut[uColIndexOut];

	if (g_bNormalizeCounts)
		{
		const FCOUNT fA = PPA.m_fOcc*wA/(wA + wB);
		const FCOUNT fB = PPB.m_fOcc*wB/(wA + wB);
		FCOUNT fTotal = 0;
		for (unsigned i = 0; i < 20; ++i)
			{
			const FCOUNT f = fA*PPA.m_fcCounts[i] + fB*PPB.m_fcCounts[i];
			PPO.m_fcCounts[i] = f;
			fTotal += f;
			}
		if (fTotal > 0)
			for (unsigned i = 0; i < 20; ++i)
				PPO.m_fcCounts[i] /= fTotal;
		}
	else
		{
		for (unsigned i = 0; i < 20; ++i)
			PPO.m_fcCounts[i] = wA*PPA.m_fcCounts[i] + wB*PPB.m_fcCounts[i];
		}
	}

void AlignTwoProfsGivenPath(const PWPath &Path,
  const ProfPos *PA, unsigned uPrefixLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uPrefixLengthB, WEIGHT wB,
  ProfPos **ptrPOut, unsigned *ptruLengthOut)
	{
#if	TRACE
	Log("AlignTwoProfsGivenPath wA=%.3g wB=%.3g Path=\n", wA, wB);
	Path.LogMe();
#endif
	assert(BTEq(wA + wB, 1.0));

	unsigned uColIndexA = 0;
	unsigned uColIndexB = 0;
	unsigned uColIndexOut = 0;
	const unsigned uEdgeCount = Path.GetEdgeCount();
	ProfPos *POut = new ProfPos[uEdgeCount];
	char cPrevType = 'M';
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		const char cType = Edge.cType;

		const unsigned uPrefixLengthA = Edge.uPrefixLengthA;
		const unsigned uPrefixLengthB = Edge.uPrefixLengthB;

#if	TRACE
		Log("\nEdge %u %c%u.%u ColA=%u ColB=%u\n",
		  uEdgeIndex,
		  Edge.cType,
		  Edge.uPrefixLengthA,
		  Edge.uPrefixLengthB,
		  uColIndexA,
		  uColIndexB);
#endif

		POut[uColIndexOut].m_bAllGaps = false;
		switch (cType)
			{
		case 'M':
			{
			assert(uPrefixLengthA > 0);
			assert(uPrefixLengthB > 0);
			SetFreqs(
			  PA, uPrefixLengthA, wA,
			  PB, uPrefixLengthB, wB,
			  POut, uColIndexOut);
			switch (cPrevType)
				{
			case 'M':
				SetGapsMM(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
			  break;
			case 'D':
				SetGapsDM(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			case 'I':
				SetGapsIM(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			default:
				Quit("Bad cPrevType");
				}
			++uColIndexA;
			++uColIndexB;
			++uColIndexOut;
			break;
			}
		case 'D':
			{
			assert(uPrefixLengthA > 0);
			SetFreqs(
			  PA, uPrefixLengthA, wA,
			  PB, uPrefixLengthB, 0,
			  POut, uColIndexOut);
			switch (cPrevType)
				{
			case 'M':
				SetGapsMD(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
			  break;
			case 'D':
				SetGapsDD(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			case 'I':
				SetGapsID(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			default:
				Quit("Bad cPrevType");
				}
			++uColIndexA;
			++uColIndexOut;
			break;
			}
		case 'I':
			{
			assert(uPrefixLengthB > 0);
			SetFreqs(
			  PA, uPrefixLengthA, 0,
			  PB, uPrefixLengthB, wB,
			  POut, uColIndexOut);
			switch (cPrevType)
				{
			case 'M':
				SetGapsMI(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
			  break;
			case 'D':
				SetGapsDI(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			case 'I':
				SetGapsII(
				  PA, uPrefixLengthA, wA,
				  PB, uPrefixLengthB, wB,
				  POut, uColIndexOut);
				break;
			default:
				Quit("Bad cPrevType");
				}
			++uColIndexB;
			++uColIndexOut;
			break;
			}
		default:
			assert(false);
			}
		cPrevType = cType;
		}
	assert(uColIndexOut == uEdgeCount);

	ProfScoresFromFreqs(POut, uEdgeCount);
	ValidateProf(POut, uEdgeCount);

	*ptrPOut = POut;
	*ptruLengthOut = uEdgeCount;

#if	TRACE
	Log("AlignTwoProfsGivenPath:\n");
	ListProfile(POut, uEdgeCount, 0);
#endif
	}
