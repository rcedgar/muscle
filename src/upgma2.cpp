#include "muscle.h"
#include "tree.h"
#include "distcalc.h"

// UPGMA clustering in O(N^2) time and space.

#define	TRACE	0

#define	MIN(x, y)	((x) < (y) ? (x) : (y))
#define	MAX(x, y)	((x) > (y) ? (x) : (y))
#define	AVG(x, y)	(((x) + (y))/2)

static unsigned g_uLeafCount;
static unsigned g_uTriangleSize;
static unsigned g_uInternalNodeCount;
static unsigned g_uInternalNodeIndex;

// Triangular distance matrix is g_Dist, which is allocated
// as a one-dimensional vector of length g_uTriangleSize.
// TriangleSubscript(i,j) maps row,column=i,j to the subscript
// into this vector.
// Row / column coordinates are a bit messy.
// Initially they are leaf indexes 0..N-1.
// But each time we create a new node (=new cluster, new subtree),
// we re-use one of the two rows that become available (the children
// of the new node). This saves memory.
// We keep track of this through the g_uNodeIndex vector.
static dist_t *g_Dist;

// Distance to nearest neighbor in row i of distance matrix.
// Subscript is distance matrix row.
static dist_t *g_MinDist;

// Nearest neighbor to row i of distance matrix.
// Subscript is distance matrix row.
static unsigned *g_uNearestNeighbor;

// Node index of row i in distance matrix.
// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
// Subscript is distance matrix row.
static unsigned *g_uNodeIndex;

// The following vectors are defined on internal nodes,
// subscripts are internal node index 0..N-2.
// For g_uLeft/Right, value is the node index 0 .. 2N-2
// because a child can be internal or leaf.
static unsigned *g_uLeft;
static unsigned *g_uRight;
static dist_t *g_Height;
static dist_t *g_LeftLength;
static dist_t *g_RightLength;

static inline unsigned TriangleSubscript(unsigned uIndex1, unsigned uIndex2)
	{
#if	DEBUG
	if (uIndex1 >= g_uLeafCount || uIndex2 >= g_uLeafCount)
		Quit("TriangleSubscript(%u,%u) %u", uIndex1, uIndex2, g_uLeafCount);
#endif
	unsigned v;
	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
	assert(v < (g_uLeafCount*(g_uLeafCount - 1))/2);
	return v;
	}

static void ListState()
	{
	Log("Dist matrix\n");
	Log("     ");
	for (unsigned i = 0; i < g_uLeafCount; ++i)
		{
		if (uInsane == g_uNodeIndex[i])
			continue;
		Log("  %5u", g_uNodeIndex[i]);
		}
	Log("\n");

	for (unsigned i = 0; i < g_uLeafCount; ++i)
		{
		if (uInsane == g_uNodeIndex[i])
			continue;
		Log("%5u  ", g_uNodeIndex[i]);
		for (unsigned j = 0; j < g_uLeafCount; ++j)
			{
			if (uInsane == g_uNodeIndex[j])
				continue;
			if (i == j)
				Log("       ");
			else
				{
				unsigned v = TriangleSubscript(i, j);
				Log("%5.2g  ", g_Dist[v]);
				}
			}
		Log("\n");
		}

	Log("\n");
	Log("    i   Node   NrNb      Dist\n");
	Log("-----  -----  -----  --------\n");
	for (unsigned i = 0; i < g_uLeafCount; ++i)
		{
		if (uInsane == g_uNodeIndex[i])
			continue;
		Log("%5u  %5u  %5u  %8.3f\n",
		  i,
		  g_uNodeIndex[i],
		  g_uNearestNeighbor[i],
		  g_MinDist[i]);
		}

	Log("\n");
	Log(" Node      L      R  Height  LLength  RLength\n");
	Log("-----  -----  -----  ------  -------  -------\n");
	for (unsigned i = 0; i <= g_uInternalNodeIndex; ++i)
		Log("%5u  %5u  %5u  %6.2g  %6.2g  %6.2g\n",
		  i,
		  g_uLeft[i],
		  g_uRight[i],
		  g_Height[i],
		  g_LeftLength[i],
		  g_RightLength[i]);
	}

void UPGMA2(const DistCalc &DC, Tree &tree, LINKAGE Linkage)
	{
	const char *DCType = DC.GetType();
	g_uLeafCount = DC.GetCount();

	g_uTriangleSize = (g_uLeafCount*(g_uLeafCount - 1))/2;
	g_uInternalNodeCount = g_uLeafCount - 1;

	g_Dist = new dist_t[g_uTriangleSize];

	g_uNodeIndex = new unsigned[g_uLeafCount];
	g_uNearestNeighbor = new unsigned[g_uLeafCount];
	g_MinDist = new dist_t[g_uLeafCount];
	unsigned *Ids = new unsigned [g_uLeafCount];
	char **Names = new char *[g_uLeafCount];

	g_uLeft = new unsigned[g_uInternalNodeCount];
	g_uRight = new unsigned[g_uInternalNodeCount];
	g_Height = new dist_t[g_uInternalNodeCount];
	g_LeftLength = new dist_t[g_uInternalNodeCount];
	g_RightLength = new dist_t[g_uInternalNodeCount];

	for (unsigned i = 0; i < g_uLeafCount; ++i)
		{
		g_MinDist[i] = BIG_DIST;
		g_uNodeIndex[i] = i;
		g_uNearestNeighbor[i] = uInsane;
		Ids[i] = DC.GetId(i);
		Names[i] = strsave(DC.GetName(i));
		}

	for (unsigned i = 0; i < g_uInternalNodeCount; ++i)
		{
		g_uLeft[i] = uInsane;
		g_uRight[i] = uInsane;
		g_LeftLength[i] = BIG_DIST;
		g_RightLength[i] = BIG_DIST;
		g_Height[i] = BIG_DIST;
		}

// Compute initial NxN triangular distance matrix.
// Store minimum distance for each full (not triangular) row.
// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
// so nothing to do when i=0.
	uint ThreadCount = GetRequestedThreadCount();
	uint Counter = 0;
	omp_lock_t Lock;
	omp_init_lock(&Lock);
#pragma omp parallel for num_threads(ThreadCount)
	for (int i = 1; i < (int) g_uLeafCount; ++i)
		{
		omp_set_lock(&Lock);
		ProgressStep(Counter, g_uLeafCount-1, "Dist mx (%s)", DCType);
		++Counter;
		omp_unset_lock(&Lock);
		dist_t *Row = g_Dist + TriangleSubscript(i, 0);
		DC.CalcDistRange(i, Row);
		for (int j = 0; j < i; ++j)
			{
			dist_t d = Row[j];
			if (d < 0)
				{
				Row[j] = 0;
				d = 0;
				}
			if (d < g_MinDist[i])
				{
				g_MinDist[i] = d;
				g_uNearestNeighbor[i] = j;
				}
			if (d < g_MinDist[j])
				{
				g_MinDist[j] = d;
				g_uNearestNeighbor[j] = i;
				}
			}
		}

#if	TRACE
	Log("Initial state:\n");
	ListState();
#endif

	for (g_uInternalNodeIndex = 0; g_uInternalNodeIndex < g_uLeafCount - 1;
	  ++g_uInternalNodeIndex)
		{
		ProgressStep(g_uInternalNodeIndex, g_uLeafCount - 1, "UPGMA2");
#if	TRACE
		Log("\n");
		Log("Internal node index %5u\n", g_uInternalNodeIndex);
		Log("-------------------------\n");
#endif

	// Find nearest neighbors
		unsigned Lmin = uInsane;
		unsigned Rmin = uInsane;
		dist_t dtMinDist = BIG_DIST;
		for (unsigned j = 0; j < g_uLeafCount; ++j)
			{
			if (uInsane == g_uNodeIndex[j])
				continue;

			dist_t d = g_MinDist[j];
			if (d < dtMinDist)
				{
				dtMinDist = d;
				Lmin = j;
				Rmin = g_uNearestNeighbor[j];
				assert(uInsane != Rmin);
				assert(uInsane != g_uNodeIndex[Rmin]);
				}
			}

		assert(Lmin != uInsane);
		assert(Rmin != uInsane);
		assert(dtMinDist != BIG_DIST);

#if	TRACE
		Log("Nearest neighbors Lmin %u[=%u] Rmin %u[=%u] dist %.3g\n",
		  Lmin,
		  g_uNodeIndex[Lmin],
		  Rmin,
		  g_uNodeIndex[Rmin],
		  dtMinDist);
#endif

	// Compute distances to new node
	// New node overwrites row currently assigned to Lmin
		dist_t dtNewMinDist = BIG_DIST;
		unsigned uNewNearestNeighbor = uInsane;
		for (unsigned j = 0; j < g_uLeafCount; ++j)
			{
			if (j == Lmin || j == Rmin)
				continue;
			if (uInsane == g_uNodeIndex[j])
				continue;

			const unsigned vL = TriangleSubscript(Lmin, j);
			const unsigned vR = TriangleSubscript(Rmin, j);
			const dist_t dL = g_Dist[vL];
			const dist_t dR = g_Dist[vR];
			dist_t dtNewDist;

			switch (Linkage)
				{
			case LINKAGE_Avg:
				dtNewDist = AVG(dL, dR);
				break;

			case LINKAGE_Min:
				dtNewDist = MIN(dL, dR);
				break;

			case LINKAGE_Max:
				dtNewDist = MAX(dL, dR);
				break;

			case LINKAGE_Biased:
				dtNewDist = g_dSUEFF*AVG(dL, dR) + (1 - g_dSUEFF)*MIN(dL, dR);
				break;

			default:
				Quit("UPGMA2: Invalid LINKAGE_%u", Linkage);
				}

		// Nasty special case.
		// If nearest neighbor of j is Lmin or Rmin, then make the new
		// node (which overwrites the row currently occupied by Lmin)
		// the nearest neighbor. This situation can occur when there are
		// equal distances in the matrix. If we don't make this fix,
		// the nearest neighbor pointer for j would become invalid.
		// (We don't need to test for == Lmin, because in that case
		// the net change needed is zero due to the change in row
		// numbering).
			if (g_uNearestNeighbor[j] == Rmin)
				g_uNearestNeighbor[j] = Lmin;

#if	TRACE
			Log("New dist to %u = (%u/%.3g + %u/%.3g)/2 = %.3g\n",
			  j, Lmin, dL, Rmin, dR, dtNewDist);
#endif
			g_Dist[vL] = dtNewDist;
			if (dtNewDist < dtNewMinDist)
				{
				dtNewMinDist = dtNewDist;
				uNewNearestNeighbor = j;
				}
			}

		assert(g_uInternalNodeIndex < g_uLeafCount - 1 || BIG_DIST != dtNewMinDist);
		assert(g_uInternalNodeIndex < g_uLeafCount - 1 || uInsane != uNewNearestNeighbor);

		const unsigned v = TriangleSubscript(Lmin, Rmin);
		const dist_t dLR = g_Dist[v];
		const dist_t dHeightNew = dLR/2;
		const unsigned uLeft = g_uNodeIndex[Lmin];
		const unsigned uRight = g_uNodeIndex[Rmin];
		const dist_t HeightLeft =
		  uLeft < g_uLeafCount ? 0 : g_Height[uLeft - g_uLeafCount];
		const dist_t HeightRight =
		  uRight < g_uLeafCount ? 0 : g_Height[uRight - g_uLeafCount];

		g_uLeft[g_uInternalNodeIndex] = uLeft;
		g_uRight[g_uInternalNodeIndex] = uRight;
		g_LeftLength[g_uInternalNodeIndex] = dHeightNew - HeightLeft;
		g_RightLength[g_uInternalNodeIndex] = dHeightNew - HeightRight;
		g_Height[g_uInternalNodeIndex] = dHeightNew;

	// Row for left child overwritten by row for new node
		g_uNodeIndex[Lmin] = g_uLeafCount + g_uInternalNodeIndex;
		g_uNearestNeighbor[Lmin] = uNewNearestNeighbor;
		g_MinDist[Lmin] = dtNewMinDist;

	// Delete row for right child
		g_uNodeIndex[Rmin] = uInsane;

#if	TRACE
		Log("\nInternalNodeIndex=%u Lmin=%u Rmin=%u\n",
		  g_uInternalNodeIndex, Lmin, Rmin);
		ListState();
#endif
		}

	unsigned uRoot = g_uLeafCount - 2;
	tree.Create(g_uLeafCount, uRoot, g_uLeft, g_uRight, g_LeftLength, g_RightLength,
	  Ids, Names);

#if	TRACE
	tree.LogMe();
#endif

	delete[] g_Dist;

	delete[] g_uNodeIndex;
	delete[] g_uNearestNeighbor;
	delete[] g_MinDist;
	delete[] g_Height;

	delete[] g_uLeft;
	delete[] g_uRight;
	delete[] g_LeftLength;
	delete[] g_RightLength;
	
	for (unsigned i = 0; i < g_uLeafCount; ++i)
		free(Names[i]);
	delete[] Names;
	delete[] Ids;
	}

class DistCalcTest : public DistCalc
	{
	virtual void CalcDistRange(unsigned i, dist_t Dist[]) const
		{
		static dist_t TestDist[5][5] =
			{
			0,  2, 14, 14, 20,
			2,  0, 14, 14, 20,
			14, 14,  0,  4, 20,
			14, 14,  4,  0, 20,
			20, 20, 20, 20,  0,
			};
		for (unsigned j = 0; j < i; ++j)
			Dist[j] = TestDist[i][j];
		}
	virtual unsigned GetCount() const
		{
		return 5;
		}
	virtual unsigned GetId(unsigned i) const
		{
		return i;
		}
	virtual const char *GetName(unsigned i) const
		{
		return "name";
		}
	};

void Test()
	{
	SetListFileName("c:\\tmp\\lobster.log", false);
	DistCalcTest DC;
	Tree tree;
	UPGMA2(DC, tree, LINKAGE_Avg);
	}
