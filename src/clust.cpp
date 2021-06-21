#include "muscle.h"
#include "clust.h"
#include "clustset.h"

#define TRACE		0

Clust::Clust()
	{
	m_Nodes = 0;
	m_uNodeCount = 0;
	m_uLeafCount = 0;
	m_uClusterCount = 0;
	m_JoinStyle = JOIN_Undefined;
	m_dDist = 0;
	m_uLeafCount = 0;
	m_ptrSet = 0;
	}

Clust::~Clust()
	{
	delete[] m_Nodes;
	delete[] m_dDist;
	delete[] m_ClusterIndexToNodeIndex;
	}

void Clust::Create(ClustSet &Set, CLUSTER Method)
	{
	m_ptrSet = &Set;

	SetLeafCount(Set.GetLeafCount());

	switch (Method)
		{
	case CLUSTER_UPGMA:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Avg;
		break;

	case CLUSTER_UPGMAMax:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Max;
		break;

	case CLUSTER_UPGMAMin:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Min;
		break;

	case CLUSTER_UPGMB:
		m_JoinStyle = JOIN_NearestNeighbor;
		m_CentroidStyle = LINKAGE_Biased;
		break;

	case CLUSTER_NeighborJoining:
		m_JoinStyle = JOIN_NeighborJoining;
		m_CentroidStyle = LINKAGE_NeighborJoining;
		break;

	default:
		Quit("Clust::Create, invalid method %d", Method);
		}

	if (m_uLeafCount <= 1)
		Quit("Clust::Create: no leaves");

	m_uNodeCount = 2*m_uLeafCount - 1;
	m_Nodes = new ClustNode[m_uNodeCount];
	m_ClusterIndexToNodeIndex = new unsigned[m_uLeafCount];

	m_ptrClusterList = 0;
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		ClustNode &Node = m_Nodes[uNodeIndex];
		Node.m_uIndex = uNodeIndex;
		if (uNodeIndex < m_uLeafCount)
			{
			Node.m_uSize = 1;
			Node.m_uLeafIndexes = new unsigned[1];
			Node.m_uLeafIndexes[0] = uNodeIndex;
			AddToClusterList(uNodeIndex);
			}
		else
			Node.m_uSize = 0;
		}

// Compute initial distance matrix between leaves
	SetProgressDesc("Build dist matrix");
	unsigned uPairIndex = 0;
	const unsigned uPairCount = (m_uLeafCount*(m_uLeafCount - 1))/2;
	for (unsigned i = 0; i < m_uLeafCount; ++i)
		for (unsigned j = 0; j < i; ++j)
			{
			const float dDist = (float) m_ptrSet->ComputeDist(*this, i, j);
			SetDist(i, j, dDist);
			if (0 == uPairIndex%10000)
				Progress(uPairIndex, uPairCount);
			++uPairIndex;
			}
	ProgressStepsDone();

// Call CreateCluster once for each internal node in the tree
	SetProgressDesc("Build guide tree");
	m_uClusterCount = m_uLeafCount;
	const unsigned uInternalNodeCount = m_uNodeCount - m_uLeafCount;
	for (unsigned uNodeIndex = m_uLeafCount; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		unsigned i = uNodeIndex + 1 - m_uLeafCount;
		Progress(i, uInternalNodeCount);
		CreateCluster();
		}
	ProgressStepsDone();
	}

void Clust::CreateCluster()
	{
	unsigned uLeftNodeIndex;
	unsigned uRightNodeIndex;
	float dLeftLength;
	float dRightLength;
	ChooseJoin(&uLeftNodeIndex, &uRightNodeIndex, &dLeftLength, &dRightLength);

	const unsigned uNewNodeIndex = m_uNodeCount - m_uClusterCount + 1;

	JoinNodes(uLeftNodeIndex, uRightNodeIndex, dLeftLength, dRightLength,
	  uNewNodeIndex);

#if	TRACE
	Log("Merge New=%u L=%u R=%u Ld=%7.2g Rd=%7.2g\n",
	  uNewNodeIndex, uLeftNodeIndex, uRightNodeIndex, dLeftLength, dRightLength);
#endif

// Compute distances to other clusters
	--m_uClusterCount;
	for (unsigned uNodeIndex = GetFirstCluster(); uNodeIndex != uInsane;
	  uNodeIndex = GetNextCluster(uNodeIndex))
		{
		if (uNodeIndex == uLeftNodeIndex || uNodeIndex == uRightNodeIndex)
			continue;

		if (uNewNodeIndex == uNodeIndex)
			continue;

		const float dDist = ComputeDist(uNewNodeIndex, uNodeIndex);
		SetDist(uNewNodeIndex, uNodeIndex, dDist);
		}

	for (unsigned uNodeIndex = GetFirstCluster(); uNodeIndex != uInsane;
	  uNodeIndex = GetNextCluster(uNodeIndex))
		{
		if (uNodeIndex == uLeftNodeIndex || uNodeIndex == uRightNodeIndex)
			continue;

		if (uNewNodeIndex == uNodeIndex)
			continue;

#if	REDLACK
		const float dMetric = ComputeMetric(uNewNodeIndex, uNodeIndex);
		InsertMetric(uNewNodeIndex, uNodeIndex, dMetric);
#endif
		}
	}

void Clust::ChooseJoin(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
  float *ptrdLeftLength, float *ptrdRightLength)
	{
	switch (m_JoinStyle)
		{
	case JOIN_NearestNeighbor:
		ChooseJoinNearestNeighbor(ptruLeftIndex, ptruRightIndex, ptrdLeftLength,
		  ptrdRightLength);
		return;
	case JOIN_NeighborJoining:
		ChooseJoinNeighborJoining(ptruLeftIndex, ptruRightIndex, ptrdLeftLength,
		  ptrdRightLength);
		return;
		}
	Quit("Clust::ChooseJoin, Invalid join style %u", m_JoinStyle);
	}

void Clust::ChooseJoinNearestNeighbor(unsigned *ptruLeftIndex,
  unsigned *ptruRightIndex, float *ptrdLeftLength, float *ptrdRightLength)
	{
	const unsigned uClusterCount = GetClusterCount();

	unsigned uMinLeftNodeIndex;
	unsigned uMinRightNodeIndex;
	GetMinMetric(&uMinLeftNodeIndex, &uMinRightNodeIndex);

	float dMinDist = GetDist(uMinLeftNodeIndex, uMinRightNodeIndex);

	const float dLeftHeight = GetHeight(uMinLeftNodeIndex);
	const float dRightHeight = GetHeight(uMinRightNodeIndex);

	*ptruLeftIndex = uMinLeftNodeIndex;
	*ptruRightIndex = uMinRightNodeIndex;
	*ptrdLeftLength = dMinDist/2 - dLeftHeight;
	*ptrdRightLength = dMinDist/2 - dRightHeight;
	}

void Clust::ChooseJoinNeighborJoining(unsigned *ptruLeftIndex,
  unsigned *ptruRightIndex, float *ptrdLeftLength, float *ptrdRightLength)
	{
	const unsigned uClusterCount = GetClusterCount();

	//unsigned uMinLeftNodeIndex = uInsane;
	//unsigned uMinRightNodeIndex = uInsane;
	//float dMinD = PLUS_INFINITY;
	//for (unsigned i = GetFirstCluster(); i != uInsane; i = GetNextCluster(i))
	//	{
	//	const float ri = Calc_r(i);
	//	for (unsigned j = GetNextCluster(i); j != uInsane; j = GetNextCluster(j))
	//		{
	//		const float rj = Calc_r(j);
	//		const float dij = GetDist(i, j);
	//		const float Dij = dij - (ri + rj);
	//		if (Dij < dMinD)
	//			{
	//			dMinD = Dij;
	//			uMinLeftNodeIndex = i;
	//			uMinRightNodeIndex = j;
	//			}
	//		}
	//	}

	unsigned uMinLeftNodeIndex;
	unsigned uMinRightNodeIndex;
	GetMinMetric(&uMinLeftNodeIndex, &uMinRightNodeIndex);

	const float dDistLR = GetDist(uMinLeftNodeIndex, uMinRightNodeIndex);
	const float rL = Calc_r(uMinLeftNodeIndex);
	const float rR = Calc_r(uMinRightNodeIndex);

	const float dLeftLength = (dDistLR + rL - rR)/2;
	const float dRightLength = (dDistLR - rL + rR)/2;

	*ptruLeftIndex = uMinLeftNodeIndex;
	*ptruRightIndex = uMinRightNodeIndex;
	*ptrdLeftLength = dLeftLength;
	*ptrdRightLength = dRightLength;
	}

void Clust::JoinNodes(unsigned uLeftIndex, unsigned uRightIndex, float dLeftLength,
  float dRightLength, unsigned uNodeIndex)
	{
	ClustNode &Parent = m_Nodes[uNodeIndex];
	ClustNode &Left = m_Nodes[uLeftIndex];
	ClustNode &Right = m_Nodes[uRightIndex];

	Left.m_dLength = dLeftLength;
	Right.m_dLength = dRightLength;

	Parent.m_ptrLeft = &Left;
	Parent.m_ptrRight = &Right;

	Left.m_ptrParent = &Parent;
	Right.m_ptrParent = &Parent;

	const unsigned uLeftSize = Left.m_uSize;
	const unsigned uRightSize = Right.m_uSize;
	const unsigned uParentSize = uLeftSize + uRightSize;
	Parent.m_uSize = uParentSize;

	assert(0 == Parent.m_uLeafIndexes);
	Parent.m_uLeafIndexes = new unsigned[uParentSize];

	const unsigned uLeftBytes = uLeftSize*sizeof(unsigned);
	const unsigned uRightBytes = uRightSize*sizeof(unsigned);
	memcpy(Parent.m_uLeafIndexes, Left.m_uLeafIndexes, uLeftBytes);
	memcpy(Parent.m_uLeafIndexes + uLeftSize, Right.m_uLeafIndexes, uRightBytes);

	DeleteFromClusterList(uLeftIndex);
	DeleteFromClusterList(uRightIndex);
	AddToClusterList(uNodeIndex);
	}

float Clust::Calc_r(unsigned uNodeIndex) const
	{
	const unsigned uClusterCount = GetClusterCount();
	if (2 == uClusterCount)
		return 0;

	float dSum = 0;
	for (unsigned i = GetFirstCluster(); i != uInsane; i = GetNextCluster(i))
		{
		if (i == uNodeIndex)
			continue;
		dSum += GetDist(uNodeIndex, i);
		}
	return dSum/(uClusterCount - 2);
	}

float Clust::ComputeDist(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	switch (m_CentroidStyle)
		{
	case LINKAGE_Avg:
		return ComputeDistAverageLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Min:
		return ComputeDistMinLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Max:
		return ComputeDistMaxLinkage(uNewNodeIndex, uNodeIndex);

	case LINKAGE_Biased:
		return ComputeDistMAFFT(uNewNodeIndex, uNodeIndex);

	case LINKAGE_NeighborJoining:
		return ComputeDistNeighborJoining(uNewNodeIndex, uNodeIndex);
		}
	Quit("Clust::ComputeDist, invalid centroid style %u", m_CentroidStyle);
	return (float) g_dNAN;
	}

float Clust::ComputeDistMinLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const float dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const float dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	return (dDistL < dDistR ? dDistL : dDistR);
	}

float Clust::ComputeDistMaxLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const float dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const float dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	return (dDistL > dDistR ? dDistL : dDistR);
	}

float Clust::ComputeDistAverageLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const float dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const float dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	return (dDistL + dDistR)/2;
	}

float Clust::ComputeDistNeighborJoining(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);
	const float dDistLR = GetDist(uLeftNodeIndex, uRightNodeIndex);
	const float dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const float dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	const float dDist = (dDistL + dDistR - dDistLR)/2;
	return dDist;
	}

// This is a mysterious variant of UPGMA reverse-engineered from MAFFT source.
float Clust::ComputeDistMAFFT(unsigned uNewNodeIndex, unsigned uNodeIndex)
	{
	const unsigned uLeftNodeIndex = GetLeftIndex(uNewNodeIndex);
	const unsigned uRightNodeIndex = GetRightIndex(uNewNodeIndex);

	const float dDistLR = GetDist(uLeftNodeIndex, uRightNodeIndex);
	const float dDistL = GetDist(uLeftNodeIndex, uNodeIndex);
	const float dDistR = GetDist(uRightNodeIndex, uNodeIndex);
	const float dMinDistLR = (dDistL < dDistR ? dDistL : dDistR);
	const float dSumDistLR = dDistL + dDistR;
	const float dDist = dMinDistLR*(1 - g_dSUEFF) + dSumDistLR*g_dSUEFF/2;
	return dDist;
	}

unsigned Clust::GetClusterCount() const
	{
	return m_uClusterCount;
	}

void Clust::LogMe() const
	{
	Log("Clust %u leaves, %u nodes, %u clusters.\n",
	  m_uLeafCount, m_uNodeCount, m_uClusterCount);

	Log("Distance matrix\n");
	const unsigned uNodeCount = GetNodeCount();
	Log("       ");
	for (unsigned i = 0; i < uNodeCount - 1; ++i)
		Log(" %7u", i);
	Log("\n");

	Log("       ");
	for (unsigned i = 0; i < uNodeCount - 1; ++i)
		Log("  ------");
	Log("\n");

	for (unsigned i = 0; i < uNodeCount - 1; ++i)
		{
		Log("%4u:  ", i);
		for (unsigned j = 0; j < i; ++j)
			Log(" %7.2g", GetDist(i, j));
		Log("\n");
		}

	Log("\n");
	Log("Node  Size  Prnt  Left  Rght   Length  Name\n");
	Log("----  ----  ----  ----  ----   ------  ----\n");
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		const ClustNode &Node = m_Nodes[uNodeIndex];
		Log("%4u  %4u", uNodeIndex, Node.m_uSize);
		if (0 != Node.m_ptrParent)
			Log("  %4u", Node.m_ptrParent->m_uIndex);
		else
			Log("      ");

		if (0 != Node.m_ptrLeft)
			Log("  %4u", Node.m_ptrLeft->m_uIndex);
		else
			Log("      ");

		if (0 != Node.m_ptrRight)
			Log("  %4u", Node.m_ptrRight->m_uIndex);
		else
			Log("      ");

		if (uNodeIndex != m_uNodeCount - 1)
			Log("  %7.3g", Node.m_dLength);
		if (IsLeaf(uNodeIndex))
			{
			const char *ptrName = GetNodeName(uNodeIndex);
			if (0 != ptrName)
				Log("  %s", ptrName);
			}
		if (GetRootNodeIndex() == uNodeIndex)
			Log("    [ROOT]");
		Log("\n");
		}
	}

const ClustNode &Clust::GetNode(unsigned uNodeIndex) const
	{
	if (uNodeIndex >= m_uNodeCount)
		Quit("ClustNode::GetNode(%u) %u", uNodeIndex, m_uNodeCount);
	return m_Nodes[uNodeIndex];
	}

bool Clust::IsLeaf(unsigned uNodeIndex) const
	{
	return uNodeIndex < m_uLeafCount;
	}

unsigned Clust::GetClusterSize(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	return Node.m_uSize;
	}

unsigned Clust::GetLeftIndex(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	if (0 == Node.m_ptrLeft)
		Quit("Clust::GetLeftIndex: leaf");
	return Node.m_ptrLeft->m_uIndex;
	}

unsigned Clust::GetRightIndex(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	if (0 == Node.m_ptrRight)
		Quit("Clust::GetRightIndex: leaf");
	return Node.m_ptrRight->m_uIndex;
	}

float Clust::GetLength(unsigned uNodeIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	return Node.m_dLength;
	}

void Clust::SetLeafCount(unsigned uLeafCount)
	{
	if (uLeafCount <= 1)
		Quit("Clust::SetLeafCount(%u)", uLeafCount);

	m_uLeafCount = uLeafCount;
	const unsigned uNodeCount = GetNodeCount();

// Triangular matrix size excluding diagonal (all zeros in our case).
	m_uTriangularMatrixSize = (uNodeCount*(uNodeCount - 1))/2;
	m_dDist = new float[m_uTriangularMatrixSize];
	}

unsigned Clust::GetLeafCount() const
	{
	return m_uLeafCount;
	}

unsigned Clust::VectorIndex(unsigned uIndex1, unsigned uIndex2) const
	{
	const unsigned uNodeCount = GetNodeCount();
	if (uIndex1 >= uNodeCount || uIndex2 >= uNodeCount)
		Quit("DistVectorIndex(%u,%u) %u", uIndex1, uIndex2, uNodeCount);
	unsigned v;
	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
	assert(v < m_uTriangularMatrixSize);
	return v;
	}

float Clust::GetDist(unsigned uIndex1, unsigned uIndex2) const
	{
	unsigned v = VectorIndex(uIndex1, uIndex2);
	return m_dDist[v];
	}

void Clust::SetDist(unsigned uIndex1, unsigned uIndex2, float dDist)
	{
	unsigned v = VectorIndex(uIndex1, uIndex2);
	m_dDist[v] = dDist;
	}

float Clust::GetHeight(unsigned uNodeIndex) const
	{
	if (IsLeaf(uNodeIndex))
		return 0;

	const unsigned uLeftIndex = GetLeftIndex(uNodeIndex);
	const unsigned uRightIndex = GetRightIndex(uNodeIndex);
	const float dLeftLength = GetLength(uLeftIndex);
	const float dRightLength = GetLength(uRightIndex);
	const float dLeftHeight = dLeftLength + GetHeight(uLeftIndex);
	const float dRightHeight = dRightLength + GetHeight(uRightIndex);
	return (dLeftHeight + dRightHeight)/2;
	}

const char *Clust::GetNodeName(unsigned uNodeIndex) const
	{
	if (!IsLeaf(uNodeIndex))
		Quit("Clust::GetNodeName, is not leaf");
	return m_ptrSet->GetLeafName(uNodeIndex);
	}

unsigned Clust::GetNodeId(unsigned uNodeIndex) const
	{
	if (uNodeIndex >= GetLeafCount())
		return 0;
	return m_ptrSet->GetLeafId(uNodeIndex);
	}

unsigned Clust::GetLeaf(unsigned uNodeIndex, unsigned uLeafIndex) const
	{
	const ClustNode &Node = GetNode(uNodeIndex);
	const unsigned uLeafCount = Node.m_uSize;
	if (uLeafIndex >= uLeafCount)
		Quit("Clust::GetLeaf, invalid index");
	const unsigned uIndex = Node.m_uLeafIndexes[uLeafIndex];
	if (uIndex >= m_uNodeCount)
		Quit("Clust::GetLeaf, index out of range");
	return uIndex;
	}

unsigned Clust::GetFirstCluster() const
	{
	if (0 == m_ptrClusterList)
		return uInsane;
	return m_ptrClusterList->m_uIndex;
	}

unsigned Clust::GetNextCluster(unsigned uIndex) const
	{
	ClustNode *ptrNode = &m_Nodes[uIndex];
	if (0 == ptrNode->m_ptrNextCluster)
		return uInsane;
	return ptrNode->m_ptrNextCluster->m_uIndex;
	}

void Clust::DeleteFromClusterList(unsigned uNodeIndex)
	{
	assert(uNodeIndex < m_uNodeCount);
	ClustNode *ptrNode = &m_Nodes[uNodeIndex];
	ClustNode *ptrPrev = ptrNode->m_ptrPrevCluster;
	ClustNode *ptrNext = ptrNode->m_ptrNextCluster;

	if (0 != ptrNext)
		ptrNext->m_ptrPrevCluster = ptrPrev;
	if (0 == ptrPrev)
		{
		assert(m_ptrClusterList == ptrNode);
		m_ptrClusterList = ptrNext;
		}
	else
		ptrPrev->m_ptrNextCluster = ptrNext;

	ptrNode->m_ptrNextCluster = 0;
	ptrNode->m_ptrPrevCluster = 0;
	}

void Clust::AddToClusterList(unsigned uNodeIndex)
	{
	assert(uNodeIndex < m_uNodeCount);
	ClustNode *ptrNode = &m_Nodes[uNodeIndex];

	if (0 != m_ptrClusterList)
		m_ptrClusterList->m_ptrPrevCluster = ptrNode;

	ptrNode->m_ptrNextCluster = m_ptrClusterList;
	ptrNode->m_ptrPrevCluster = 0;

	m_ptrClusterList = ptrNode;
	}

float Clust::ComputeMetric(unsigned uIndex1, unsigned uIndex2) const
	{
	switch (m_JoinStyle)
		{
	case JOIN_NearestNeighbor:
		return ComputeMetricNearestNeighbor(uIndex1, uIndex2);

	case JOIN_NeighborJoining:
		return ComputeMetricNeighborJoining(uIndex1, uIndex2);
		}
	Quit("Clust::ComputeMetric");
	return 0;
	}

float Clust::ComputeMetricNeighborJoining(unsigned i, unsigned j) const
	{
	float ri = Calc_r(i);
	float rj = Calc_r(j);
	float dij = GetDist(i, j);
	float dMetric = dij - (ri + rj);
	return (float) dMetric;
	}

float Clust::ComputeMetricNearestNeighbor(unsigned i, unsigned j) const
	{
	return (float) GetDist(i, j);
	}

float Clust::GetMinMetricBruteForce(unsigned *ptruIndex1, unsigned *ptruIndex2) const
	{
	unsigned uMinLeftNodeIndex = uInsane;
	unsigned uMinRightNodeIndex = uInsane;
	float dMinMetric = PLUS_INFINITY;
	for (unsigned uLeftNodeIndex = GetFirstCluster(); uLeftNodeIndex != uInsane;
	  uLeftNodeIndex = GetNextCluster(uLeftNodeIndex))
		{
		for (unsigned uRightNodeIndex = GetNextCluster(uLeftNodeIndex);
		  uRightNodeIndex != uInsane;
		  uRightNodeIndex = GetNextCluster(uRightNodeIndex))
			{
			float dMetric = ComputeMetric(uLeftNodeIndex, uRightNodeIndex);
			if (dMetric < dMinMetric)
				{
				dMinMetric = dMetric;
				uMinLeftNodeIndex = uLeftNodeIndex;
				uMinRightNodeIndex = uRightNodeIndex;
				}
			}
		}
	*ptruIndex1 = uMinLeftNodeIndex;
	*ptruIndex2 = uMinRightNodeIndex;
	return dMinMetric;
	}

float Clust::GetMinMetric(unsigned *ptruIndex1, unsigned *ptruIndex2) const
	{
	return GetMinMetricBruteForce(ptruIndex1, ptruIndex2);
	}
