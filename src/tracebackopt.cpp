#include "muscle.h"
#include "pwpath.h"

void TraceBackToPath(int **TraceBack, unsigned uLengthA,
  unsigned uLengthB, PWPath &Path)
	{
	Path.Clear();

	PWEdge Edge;
	Edge.uPrefixLengthA = uLengthA;
	Edge.uPrefixLengthB = uLengthB;

	for (;;)
		{
		if (0 == Edge.uPrefixLengthA && 0 == Edge.uPrefixLengthB)
			break;

		int iDelta = TraceBack[Edge.uPrefixLengthA][Edge.uPrefixLengthB];
#if	TRACE
		Log("TraceBack[%u][%u] = %d\n",
		  Edge.uPrefixLengthA, Edge.uPrefixLengthB, iDelta);
#endif
		if (0 == iDelta)
			{
			assert(Edge.uPrefixLengthA > 0);
			assert(Edge.uPrefixLengthB > 0);

			Edge.cType = 'M';
			Path.PrependEdge(Edge);
			--(Edge.uPrefixLengthA);
			--(Edge.uPrefixLengthB);
			continue;
			}
		else if (iDelta > 0)
			{
			Edge.cType = 'D';
			while (iDelta-- > 0)
				{
				assert(Edge.uPrefixLengthA > 0);

				Path.PrependEdge(Edge);
				--(Edge.uPrefixLengthA);
				}
			}
		else if (iDelta < 0)
			{
			Edge.cType = 'I';
			while (iDelta++ < 0)
				{
				assert(Edge.uPrefixLengthB > 0);

				Path.PrependEdge(Edge);
				--(Edge.uPrefixLengthB);
				}
			}

		if (0 == Edge.uPrefixLengthA && 0 == Edge.uPrefixLengthB)
			break;

		assert(Edge.uPrefixLengthA > 0);
		assert(Edge.uPrefixLengthB > 0);

		Edge.cType = 'M';
		Path.PrependEdge(Edge);
		--(Edge.uPrefixLengthA);
		--(Edge.uPrefixLengthB);
		}

#if	TRACE
	Log("TraceBackToPath ");
	Path.LogMe();
#endif
	}
