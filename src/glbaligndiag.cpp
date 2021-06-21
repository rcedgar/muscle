#include "muscle.h"
#include "dpreglist.h"
#include "diaglist.h"
#include "pwpath.h"
#include "profile.h"
#include "timing.h"

#define TRACE		0
#define TRACE_PATH	0
#define LIST_DIAGS	0

static double g_dDPAreaWithoutDiags = 0.0;
static double g_dDPAreaWithDiags = 0.0;

static void OffsetPath(PWPath &Path, unsigned uOffsetA, unsigned uOffsetB)
	{
	const unsigned uEdgeCount = Path.GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);

	// Nasty hack -- poke new values back into path, circumventing class
		PWEdge &NonConstEdge = (PWEdge &) Edge;
		NonConstEdge.uPrefixLengthA += uOffsetA;
		NonConstEdge.uPrefixLengthB += uOffsetB;
		}
	}

static void DiagToPath(const Diag &d, PWPath &Path)
	{
	Path.Clear();
	const unsigned uLength = d.m_uLength;
	for (unsigned i = 0; i < uLength; ++i)
		{
		PWEdge Edge;
		Edge.cType = 'M';
		Edge.uPrefixLengthA = d.m_uStartPosA + i + 1;
		Edge.uPrefixLengthB = d.m_uStartPosB + i + 1;
		Path.AppendEdge(Edge);
		}
	}

static void AppendRegPath(PWPath &Path, const PWPath &RegPath)
	{
	const unsigned uRegEdgeCount = RegPath.GetEdgeCount();
	for (unsigned uRegEdgeIndex = 0; uRegEdgeIndex < uRegEdgeCount; ++uRegEdgeIndex)
		{
		const PWEdge &RegEdge = RegPath.GetEdge(uRegEdgeIndex);
		Path.AppendEdge(RegEdge);
		}
	}

SCORE GlobalAlignDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
#if	LIST_DIAGS
	TICKS t1 = GetClockTicks();
#endif

	DiagList DL;

	if (ALPHA_Amino == g_Alpha)
		FindDiags(PA, uLengthA, PB, uLengthB, DL);
	else if (ALPHA_DNA == g_Alpha || ALPHA_RNA == g_Alpha)
		FindDiagsNuc(PA, uLengthA, PB, uLengthB, DL);
	else
		Quit("GlobalAlignDiags: bad alpha");

#if	TRACE
	Log("GlobalAlignDiags, diag list:\n");
	DL.LogMe();
#endif

	DL.Sort();
	DL.DeleteIncompatible();

#if	TRACE
	Log("After DeleteIncompatible:\n");
	DL.LogMe();
#endif

	MergeDiags(DL);

#if	TRACE
	Log("After MergeDiags:\n");
	DL.LogMe();
#endif

	DPRegionList RL;
	DiagListToDPRegionList(DL, RL, uLengthA, uLengthB);

#if	TRACE
	Log("RegionList:\n");
	RL.LogMe();
#endif

#if	LIST_DIAGS
	{
	TICKS t2 = GetClockTicks();
	unsigned uArea = RL.GetDPArea();
	Log("ticks=%ld\n", (long) (t2 - t1));
	Log("area=%u\n", uArea);
	}
#endif

	g_dDPAreaWithoutDiags += uLengthA*uLengthB;

	double dDPAreaWithDiags = 0.0;
	const unsigned uRegionCount = RL.GetCount();
	for (unsigned uRegionIndex = 0; uRegionIndex < uRegionCount; ++uRegionIndex)
		{
		const DPRegion &r = RL.Get(uRegionIndex);

		PWPath RegPath;
		if (DPREGIONTYPE_Diag == r.m_Type)
			{
			DiagToPath(r.m_Diag, RegPath);
#if	TRACE_PATH
			Log("DiagToPath, path=\n");
			RegPath.LogMe();
#endif
			}
		else if (DPREGIONTYPE_Rect == r.m_Type)
			{
			const unsigned uRegStartPosA = r.m_Rect.m_uStartPosA;
			const unsigned uRegStartPosB = r.m_Rect.m_uStartPosB;
			const unsigned uRegLengthA = r.m_Rect.m_uLengthA;
			const unsigned uRegLengthB = r.m_Rect.m_uLengthB;
			const ProfPos *RegPA = PA + uRegStartPosA;
			const ProfPos *RegPB = PB + uRegStartPosB;

			dDPAreaWithDiags += uRegLengthA*uRegLengthB;
			GlobalAlignNoDiags(RegPA, uRegLengthA, RegPB, uRegLengthB, RegPath);
#if	TRACE_PATH
			Log("GlobalAlignNoDiags RegPath=\n");
			RegPath.LogMe();
#endif
			OffsetPath(RegPath, uRegStartPosA, uRegStartPosB);
#if	TRACE_PATH
			Log("After offset path, RegPath=\n");
			RegPath.LogMe();
#endif
			}
		else
			Quit("GlobalAlignDiags, Invalid region type %u", r.m_Type);

		AppendRegPath(Path, RegPath);
#if	TRACE_PATH
		Log("After AppendPath, path=");
		Path.LogMe();
#endif
		}

#if	TRACE
	{
	double dDPAreaWithoutDiags = uLengthA*uLengthB;
	Log("DP area with diags %.3g without %.3g pct saved %.3g %%\n",
	  dDPAreaWithDiags, dDPAreaWithoutDiags, (1.0 - dDPAreaWithDiags/dDPAreaWithoutDiags)*100.0);
	}
#endif
	g_dDPAreaWithDiags += dDPAreaWithDiags;
	return 0;
	}

//void ListDiagSavings()
//	{
//	if (!g_bVerbose || !g_bDiags)
//		return;
//	double dAreaSaved = g_dDPAreaWithoutDiags - g_dDPAreaWithDiags;
//	double dPct = dAreaSaved*100.0/g_dDPAreaWithoutDiags;
//	Log("DP area saved by diagonals %-4.1f%%\n", dPct);
//	}
