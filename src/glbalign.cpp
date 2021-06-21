#include "muscle.h"
#include "pwpath.h"
#include "timing.h"
#include "textfile.h"
#include "msa.h"
#include "profile.h"

#if	!VER_3_52

#define COMPARE_SIMPLE	0

#if	TIMING
TICKS g_ticksDP = 0;
#endif

#if	1
extern bool g_bKeepSimpleDP;
SCORE NWSmall(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE NWDASmall(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE NWDASimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE NWDASimple2(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE GlobalAlignSimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);

SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	return GlobalAlign(PA, uLengthA, PB, uLengthB, Path);
	}

#if	COMPARE_SIMPLE

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	g_bKeepSimpleDP = true;
	PWPath SimplePath;
	GlobalAlignSimple(PA, uLengthA, PB, uLengthB, SimplePath);

	SCORE Score = NWSmall(PA, uLengthA, PB, uLengthB, Path);

	if (!Path.Equal(SimplePath))
		{
		Log("Simple:\n");
		SimplePath.LogMe();
		Log("Small:\n");
		Path.LogMe();
		Quit("Paths differ");
		}

#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksDP += (t2 - t1);
#endif
	return Score;
	}

#else // COMPARE_SIMPLE

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	SCORE Score = NWSmall(PA, uLengthA, PB, uLengthB, Path);
#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksDP += (t2 - t1);
#endif
	return Score;
	}

#endif

#else // 1

static void AllInserts(PWPath &Path, unsigned uLengthB)
	{
	Path.Clear();
	PWEdge Edge;
	Edge.cType = 'I';
	Edge.uPrefixLengthA = 0;
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB <= uLengthB; ++uPrefixLengthB)
		{
		Edge.uPrefixLengthB = uPrefixLengthB;
		Path.AppendEdge(Edge);
		}
	}

static void AllDeletes(PWPath &Path, unsigned uLengthA)
	{
	Path.Clear();
	PWEdge Edge;
	Edge.cType = 'D';
	Edge.uPrefixLengthB = 0;
	for (unsigned uPrefixLengthA = 1; uPrefixLengthA <= uLengthA; ++uPrefixLengthA)
		{
		Edge.uPrefixLengthA = uPrefixLengthA;
		Path.AppendEdge(Edge);
		}
	}

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	if (0 == uLengthA)
		{
		AllInserts(Path, uLengthB);
		return 0;
		}
	else if (0 == uLengthB)
		{
		AllDeletes(Path, uLengthA);
		return 0;
		}

	SCORE Score = 0;
	if (g_bDiags)
		Score = GlobalAlignDiags(PA, uLengthA, PB, uLengthB, Path);
	else
		Score = GlobalAlignNoDiags(PA, uLengthA, PB, uLengthB, Path);
#if	TIMING
	TICKS t2 = GetClockTicks();
	g_ticksDP += (t2 - t1);
#endif
	return Score;
	}

SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (g_bDimer)
		return GlobalAlignDimer(PA, uLengthA, PB, uLengthB, Path);

	switch (g_PPScore)
		{
	case PPSCORE_LE:
		return GlobalAlignLE(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SP:
	case PPSCORE_SV:
		return GlobalAlignSP(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SPN:
		return GlobalAlignSPN(PA, uLengthA, PB, uLengthB, Path);
		}

	Quit("Invalid PP score (GlobalAlignNoDiags)");
	return 0;
	}

#endif

#endif	// !VER_3_52
