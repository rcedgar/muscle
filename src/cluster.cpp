#include "muscle.h"
#include "cluster.h"
#include "distfunc.h"

static inline float Min(float d1, float d2)
	{
	return d1 < d2 ? d1 : d2;
	}

static inline float Max(float d1, float d2)
	{
	return d1 > d2 ? d1 : d2;
	}

static inline float Mean(float d1, float d2)
	{
	return (float) ((d1 + d2)/2.0);
	}

#if	_DEBUG
void ClusterTree::Validate(unsigned uNodeCount)
	{
	unsigned n;
	ClusterNode *pNode;
	unsigned uDisjointListCount = 0;
	for (pNode = m_ptrDisjoints; pNode; pNode = pNode->GetNextDisjoint())
		{
		ClusterNode *pPrev = pNode->GetPrevDisjoint();
		ClusterNode *pNext = pNode->GetNextDisjoint();
		if (0 != pPrev)
			{
			if (pPrev->GetNextDisjoint() != pNode)
				{
				Log("Prev->This mismatch, prev=\n");
				pPrev->LogMe();
				Log("This=\n");
				pNode->LogMe();
				Quit("ClusterTree::Validate()");
				}
			}
		else
			{
			if (pNode != m_ptrDisjoints)
				{
				Log("[%u]->prev = 0 but != m_ptrDisjoints=%d\n",
				  pNode->GetIndex(),
				  m_ptrDisjoints ? m_ptrDisjoints->GetIndex() : 0xffffffff);
				pNode->LogMe();
				Quit("ClusterTree::Validate()");
				}
			}
		if (0 != pNext)
			{
			if (pNext->GetPrevDisjoint() != pNode)
				{
				Log("Next->This mismatch, next=\n");
				pNext->LogMe();
				Log("This=\n");
				pNode->LogMe();
				Quit("ClusterTree::Validate()");
				}
			}
		++uDisjointListCount;
		if (uDisjointListCount > m_uNodeCount)
			Quit("Loop in disjoint list");
		}

	unsigned uParentlessNodeCount = 0;
	for (n = 0; n < uNodeCount; ++n)
		if (0 == m_Nodes[n].GetParent())
			++uParentlessNodeCount;
	
	if (uDisjointListCount != uParentlessNodeCount)
		Quit("Disjoints = %u Parentless = %u\n", uDisjointListCount,
		  uParentlessNodeCount);
	}
#else	// !_DEBUG
#define	Validate(uNodeCount)	// empty
#endif

void ClusterNode::LogMe() const
	{
	unsigned uClusterSize = GetClusterSize();
	Log("[%02u] w=%5.3f  CW=%5.3f  LBW=%5.3f  RBW=%5.3f  LWT=%5.3f  RWT=%5.3f  L=%02d  R=%02d  P=%02d  NxDj=%02d  PvDj=%02d  Sz=%02d  {",
		m_uIndex,
		m_dWeight,
		GetClusterWeight(),
		GetLeftBranchWeight(),
		GetRightBranchWeight(),
		GetLeftWeight(),
		GetRightWeight(),
		m_ptrLeft ? m_ptrLeft->GetIndex() : 0xffffffff,
		m_ptrRight ? m_ptrRight->GetIndex() : 0xffffffff,
		m_ptrParent ? m_ptrParent->GetIndex() : 0xffffffff,
		m_ptrNextDisjoint ? m_ptrNextDisjoint->GetIndex() : 0xffffffff,
		m_ptrPrevDisjoint ? m_ptrPrevDisjoint->GetIndex() : 0xffffffff,
		uClusterSize);
	for (unsigned i = 0; i < uClusterSize; ++i)
		Log(" %u", GetClusterLeaf(i)->GetIndex());
	Log(" }\n");
	}

// How many leaves in the sub-tree under this node?
unsigned ClusterNode::GetClusterSize() const
	{
	unsigned uLeafCount = 0;

	if (0 == m_ptrLeft && 0 == m_ptrRight)
		return 1;

	if (0 != m_ptrLeft)
		uLeafCount += m_ptrLeft->GetClusterSize();
	if (0 != m_ptrRight)
		uLeafCount += m_ptrRight->GetClusterSize();
	assert(uLeafCount > 0);
	return uLeafCount;
	}

double ClusterNode::GetClusterWeight() const
	{
	double dWeight = 0.0;
	if (0 != m_ptrLeft)
		dWeight += m_ptrLeft->GetClusterWeight();
	if (0 != m_ptrRight)
		dWeight += m_ptrRight->GetClusterWeight();
	return dWeight + GetWeight();
	}

double ClusterNode::GetLeftBranchWeight() const
	{
	const ClusterNode *ptrLeft = GetLeft();
	if (0 == ptrLeft)
		return 0.0;

	return GetWeight() - ptrLeft->GetWeight();
	}

double ClusterNode::GetRightBranchWeight() const
	{
	const ClusterNode *ptrRight = GetRight();
	if (0 == ptrRight)
		return 0.0;

	return GetWeight() - ptrRight->GetWeight();
	}

double ClusterNode::GetRightWeight() const
	{
	const ClusterNode *ptrRight = GetRight();
	if (0 == ptrRight)
		return 0.0;
	return ptrRight->GetClusterWeight() + GetWeight();
	}

double ClusterNode::GetLeftWeight() const
	{
	const ClusterNode *ptrLeft = GetLeft();
	if (0 == ptrLeft)
		return 0.0;
	return ptrLeft->GetClusterWeight() + GetWeight();
	}

// Return n'th leaf in the sub-tree under this node.
const ClusterNode *ClusterNode::GetClusterLeaf(unsigned uLeafIndex) const
	{
	if (0 != m_ptrLeft)
		{
		if (0 == m_ptrRight)
			return this;

		unsigned uLeftLeafCount = m_ptrLeft->GetClusterSize();

		if (uLeafIndex < uLeftLeafCount)
			return m_ptrLeft->GetClusterLeaf(uLeafIndex);

		assert(uLeafIndex >= uLeftLeafCount);
		return m_ptrRight->GetClusterLeaf(uLeafIndex - uLeftLeafCount);
		}
	if (0 == m_ptrRight)
		return this;
	return m_ptrRight->GetClusterLeaf(uLeafIndex);
	}

void ClusterTree::DeleteFromDisjoints(ClusterNode *ptrNode)
	{
	ClusterNode *ptrPrev = ptrNode->GetPrevDisjoint();
	ClusterNode *ptrNext = ptrNode->GetNextDisjoint();

	if (0 != ptrPrev)
		ptrPrev->SetNextDisjoint(ptrNext);
	else
		m_ptrDisjoints = ptrNext;

	if (0 != ptrNext)
		ptrNext->SetPrevDisjoint(ptrPrev);

#if	_DEBUG
// not algorithmically necessary, but improves clarity
// and supports Validate().
	ptrNode->SetPrevDisjoint(0);
	ptrNode->SetNextDisjoint(0);
#endif
	}

void ClusterTree::AddToDisjoints(ClusterNode *ptrNode)
	{
	ptrNode->SetNextDisjoint(m_ptrDisjoints);
	ptrNode->SetPrevDisjoint(0);
	if (0 != m_ptrDisjoints)
		m_ptrDisjoints->SetPrevDisjoint(ptrNode);
	m_ptrDisjoints = ptrNode;
	}

ClusterTree::ClusterTree()
	{
	m_ptrDisjoints = 0;
	m_Nodes = 0;
	m_uNodeCount = 0;
	}

ClusterTree::~ClusterTree()
	{
	delete[] m_Nodes;
	}

void ClusterTree::LogMe() const
	{
	Log("Disjoints=%d\n", m_ptrDisjoints ? m_ptrDisjoints->GetIndex() : 0xffffffff);
	for (unsigned i = 0; i < m_uNodeCount; ++i)
		{
		m_Nodes[i].LogMe();
		}
	}

ClusterNode *ClusterTree::GetRoot() const
	{
	return &m_Nodes[m_uNodeCount - 1];
	}

// This is the UPGMA algorithm as described in Durbin et al. p166.
void ClusterTree::Create(const DistFunc &Dist)
	{
	unsigned i;
	m_uLeafCount = Dist.GetCount();
	m_uNodeCount = 2*m_uLeafCount - 1;

	delete[] m_Nodes;
	m_Nodes = new ClusterNode[m_uNodeCount];

	for (i = 0; i < m_uNodeCount; ++i)
		m_Nodes[i].SetIndex(i);

	for (i = 0; i < m_uLeafCount - 1; ++i)
		m_Nodes[i].SetNextDisjoint(&m_Nodes[i+1]);

	for (i = 1; i < m_uLeafCount; ++i)
		m_Nodes[i].SetPrevDisjoint(&m_Nodes[i-1]);
	
	m_ptrDisjoints = &m_Nodes[0];

//	Log("Initial state\n");
//	LogMe();
//	Log("\n");

	DistFunc ClusterDist;
	ClusterDist.SetCount(m_uNodeCount);
	double dMaxDist = 0.0;
	for (i = 0; i < m_uLeafCount; ++i)
		for (unsigned j = 0; j < m_uLeafCount; ++j)
			{
			float dDist = Dist.GetDist(i, j);
			ClusterDist.SetDist(i, j, dDist);
			}

	Validate(m_uLeafCount);

// Iteration. N-1 joins needed to create a binary tree from N leaves.
	for (unsigned uJoinIndex = m_uLeafCount; uJoinIndex < m_uNodeCount;
	  ++uJoinIndex)
		{
	// Find closest pair of clusters
		unsigned uIndexClosest1;
		unsigned uIndexClosest2;
		bool bFound = false;
		double dDistClosest = 9e99;
		for (ClusterNode *ptrNode1 = m_ptrDisjoints; ptrNode1;
		  ptrNode1 = ptrNode1->GetNextDisjoint())
			{
			for (ClusterNode *ptrNode2 = ptrNode1->GetNextDisjoint(); ptrNode2;
			  ptrNode2 = ptrNode2->GetNextDisjoint())
				{
				unsigned i1 = ptrNode1->GetIndex();
				unsigned i2 = ptrNode2->GetIndex();
				double dDist = ClusterDist.GetDist(i1, i2);
				if (dDist < dDistClosest)
					{
					bFound = true;
					dDistClosest = dDist;
					uIndexClosest1 = i1;
					uIndexClosest2 = i2;
					}
				}
			}
		assert(bFound);

		ClusterNode &Join = m_Nodes[uJoinIndex];
		ClusterNode &Child1 = m_Nodes[uIndexClosest1];
		ClusterNode &Child2 = m_Nodes[uIndexClosest2];

		Join.SetLeft(&Child1);
		Join.SetRight(&Child2);
		Join.SetWeight(dDistClosest);

		Child1.SetParent(&Join);
		Child2.SetParent(&Join);

		DeleteFromDisjoints(&Child1);
		DeleteFromDisjoints(&Child2);
		AddToDisjoints(&Join);

//		Log("After join %d %d\n", uIndexClosest1, uIndexClosest2);
//		LogMe();

	// Calculate distance of every remaining disjoint cluster to the
	// new cluster created by the join
		for (ClusterNode *ptrNode = m_ptrDisjoints; ptrNode;
		  ptrNode = ptrNode->GetNextDisjoint())
			{
			unsigned uNodeIndex = ptrNode->GetIndex();
			float dDist1 = ClusterDist.GetDist(uNodeIndex, uIndexClosest1);
			float dDist2 = ClusterDist.GetDist(uNodeIndex, uIndexClosest2);
			float dDist = Min(dDist1, dDist2);
			ClusterDist.SetDist(uJoinIndex, uNodeIndex, dDist);
			}
		Validate(uJoinIndex+1);
		}
	GetRoot()->GetClusterWeight();
//	LogMe();
	}
