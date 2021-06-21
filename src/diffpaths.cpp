#include "muscle.h"
#include "pwpath.h"

#define TRACE	0

void DiffPaths(const PWPath &p1, const PWPath &p2, unsigned Edges1[],
  unsigned *ptruDiffCount1, unsigned Edges2[], unsigned *ptruDiffCount2)
	{
#if	TRACE
	Log("DiffPaths\n");
	Log("p1=");
	p1.LogMe();
	Log("p2=");
	p2.LogMe();
#endif
	const unsigned uEdgeCount1 = p1.GetEdgeCount();
	const unsigned uEdgeCount2 = p2.GetEdgeCount();

	unsigned uDiffCount1 = 0;
	unsigned uDiffCount2 = 0;
	unsigned uEdgeIndex1 = 0;
	unsigned uEdgeIndex2 = 0;
	const PWEdge *Edge1 = &p1.GetEdge(uEdgeIndex1);
	const PWEdge *Edge2 = &p2.GetEdge(uEdgeIndex2);
	for (;;)
		{
		unsigned uEdgeIndexTop1 = uEdgeIndex1;
		unsigned uEdgeIndexTop2 = uEdgeIndex2;
		Edge1 = &p1.GetEdge(uEdgeIndex1);
		Edge2 = &p2.GetEdge(uEdgeIndex2);
#if	TRACE
		Log("e1[%u] PLA%u PLB%u %c, e2[%u] PLA%u PLB %u %c  DC1=%u DC2=%u\n",
		  uEdgeIndex1, Edge1->uPrefixLengthA, Edge1->uPrefixLengthB, Edge1->cType,
		  uEdgeIndex2, Edge2->uPrefixLengthA, Edge2->uPrefixLengthB, Edge2->cType,
		  uDiffCount1, uDiffCount2);
#endif
		if (Edge1->uPrefixLengthA == Edge2->uPrefixLengthA &&
		  Edge1->uPrefixLengthB == Edge2->uPrefixLengthB)
			{
			if (!Edge1->Equal(*Edge2))
				{
				Edges1[uDiffCount1++] = uEdgeIndex1;
				Edges2[uDiffCount2++] = uEdgeIndex2;
				}
			++uEdgeIndex1;
			++uEdgeIndex2;
			}

		else if (Edge2->uPrefixLengthA < Edge1->uPrefixLengthA ||
		  Edge2->uPrefixLengthB < Edge1->uPrefixLengthB)
			Edges2[uDiffCount2++] = uEdgeIndex2++;

		else if (Edge1->uPrefixLengthA < Edge2->uPrefixLengthA ||
		  Edge1->uPrefixLengthB < Edge2->uPrefixLengthB)
			Edges1[uDiffCount1++] = uEdgeIndex1++;

		if (uEdgeCount1 == uEdgeIndex1)
			{
			while (uEdgeIndex2 < uEdgeCount2)
				Edges2[uDiffCount2++] = uEdgeIndex2++;
			goto Done;
			}
		if (uEdgeCount2 == uEdgeIndex2)
			{
			while (uEdgeIndex1 < uEdgeCount1)
				Edges1[uDiffCount1++] = uEdgeIndex1++;
			goto Done;
			}
		if (uEdgeIndex1 == uEdgeIndexTop1 && uEdgeIndex2 == uEdgeIndexTop2)
			Quit("DiffPaths stuck");
		}
Done:;
#if	TRACE
	Log("DiffCount1=%u (%u %u)\n", uDiffCount1, uEdgeCount1, uEdgeCount2);
	Log("Diffs1=");
	for (unsigned i = 0; i < uDiffCount1; ++i)
		{
		const PWEdge e = p1.GetEdge(Edges1[i]);
		Log(" %u=%c%u.%u", Edges1[i], e.cType, e.uPrefixLengthA, e.uPrefixLengthB); 
		}
	Log("\n");
	Log("DiffCount2=%u\n", uDiffCount2);
	Log("Diffs2=");
	for (unsigned i = 0; i < uDiffCount2; ++i)
		{
		const PWEdge e = p2.GetEdge(Edges2[i]);
		Log(" %u=%c%u.%u", Edges2[i], e.cType, e.uPrefixLengthA, e.uPrefixLengthB); 
		}
	Log("\n");
#endif
	*ptruDiffCount1 = uDiffCount1;
	*ptruDiffCount2 = uDiffCount2;
	}

void TestDiffPaths()
	{
	PWPath p1;
	PWPath p2;

	p1.AppendEdge('M', 1, 1);
	p1.AppendEdge('M', 2, 2);
	p1.AppendEdge('M', 3, 3);

	p2.AppendEdge('M', 1, 1);
	p2.AppendEdge('D', 2, 1);
	p2.AppendEdge('I', 2, 2);
	p2.AppendEdge('M', 3, 3);

	unsigned Edges1[64];
	unsigned Edges2[64];
	unsigned uDiffCount1;
	unsigned uDiffCount2;
	DiffPaths(p1, p2, Edges1, &uDiffCount1, Edges2, &uDiffCount2);
	}
