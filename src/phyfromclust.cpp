#include "muscle.h"
#include "tree.h"
#include "clust.h"

void Tree::InitCache(unsigned uCacheCount)
	{
	m_uCacheCount = uCacheCount;

	m_uNeighbor1 = new unsigned[m_uCacheCount];
	m_uNeighbor2 = new unsigned[m_uCacheCount];
	m_uNeighbor3 = new unsigned[m_uCacheCount];

	m_Ids = new unsigned[m_uCacheCount];

	m_dEdgeLength1 = new double[m_uCacheCount];
	m_dEdgeLength2 = new double[m_uCacheCount];
	m_dEdgeLength3 = new double[m_uCacheCount];
	m_dHeight = new double[m_uCacheCount];

	m_bHasEdgeLength1 = new bool[m_uCacheCount];
	m_bHasEdgeLength2 = new bool[m_uCacheCount];
	m_bHasEdgeLength3 = new bool[m_uCacheCount];
	m_bHasHeight = new bool[m_uCacheCount];

	m_ptrName = new char *[m_uCacheCount];

	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		m_uNeighbor1[uNodeIndex] = NULL_NEIGHBOR;
		m_uNeighbor2[uNodeIndex] = NULL_NEIGHBOR;
		m_uNeighbor3[uNodeIndex] = NULL_NEIGHBOR;
		m_bHasEdgeLength1[uNodeIndex] = false;
		m_bHasEdgeLength2[uNodeIndex] = false;
		m_bHasEdgeLength3[uNodeIndex] = false;
		m_bHasHeight[uNodeIndex] = false;
		m_dEdgeLength1[uNodeIndex] = dInsane;
		m_dEdgeLength2[uNodeIndex] = dInsane;
		m_dEdgeLength3[uNodeIndex] = dInsane;
		m_dHeight[uNodeIndex] = dInsane;
		m_ptrName[uNodeIndex] = 0;
		m_Ids[uNodeIndex] = uInsane;
		}
	}

void Tree::FromClust(Clust &C)
	{
	Clear();

	m_uNodeCount = C.GetNodeCount();
	InitCache(m_uNodeCount);

// Cluster is always rooted. An unrooted cluster
// is represented by a pseudo-root, which we fix later.
	m_bRooted = true;
	const unsigned uRoot = C.GetRootNodeIndex();
	m_uRootNodeIndex = uRoot;
	m_uNeighbor1[uRoot] = NULL_NEIGHBOR;
	m_bHasEdgeLength1[uRoot] = false;

	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		if (C.IsLeaf(uNodeIndex))
			{
			const char *ptrName = C.GetNodeName(uNodeIndex);
			m_ptrName[uNodeIndex] = strsave(ptrName);
			m_Ids[uNodeIndex] = C.GetNodeId(uNodeIndex);
			continue;
			}

		const unsigned uLeft = C.GetLeftIndex(uNodeIndex);
		const unsigned uRight = C.GetRightIndex(uNodeIndex);

		const double dLeftLength = C.GetLength(uLeft);
		const double dRightLength = C.GetLength(uRight);

		m_uNeighbor2[uNodeIndex] = uLeft;
		m_uNeighbor3[uNodeIndex] = uRight;

		m_dEdgeLength1[uLeft] = dLeftLength;
		m_dEdgeLength1[uRight] = dRightLength;

		m_uNeighbor1[uLeft] = uNodeIndex;
		m_uNeighbor1[uRight] = uNodeIndex;

		m_bHasEdgeLength1[uLeft] = true;
		m_bHasEdgeLength1[uRight] = true;

		m_dEdgeLength2[uNodeIndex] = dLeftLength;
		m_dEdgeLength3[uNodeIndex] = dRightLength;

		m_bHasEdgeLength2[uNodeIndex] = true;
		m_bHasEdgeLength3[uNodeIndex] = true;
		}
	Validate();
	}
