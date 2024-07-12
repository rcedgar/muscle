#include "muscle.h"
#include "tree.h"
#include <math.h>

#define TRACE 0

/***
Node has 0 to 3 neighbors:
	0 neighbors:	singleton root
	1 neighbor:		leaf, neighbor is parent
	2 neigbors:		non-singleton root
	3 neighbors:	internal node (other than root)

Minimal rooted tree is single node.
Minimal unrooted tree is single edge.

Leaf node always has nulls in neighbors 2 and 3, neighbor 1 is parent.
When tree is rooted, neighbor 1=parent, 2=left, 3=right.

Tree2 from Newick
=================
    Nbr1	Nbr2	Nbr3
:-----------------------------------------------------------:
:	Nbr1	Nbr2	Nbr3	:	Non-leaf, unrooted			:
:	Parent	Left	Right	:	Internal, rooted			:
:	Parent	*		*		:	Leaf, rooted or unrooted	:
:	*		Left	Right	:	Root						:
:-----------------------------------------------------------:
***/

void Tree::InitCache(uint uCacheCount)
	{
	m_uCacheCount = uCacheCount;

	m_uNeighbor1 = myalloc(uint, m_uCacheCount);
	m_uNeighbor2 = myalloc(uint, m_uCacheCount);
	m_uNeighbor3 = myalloc(uint, m_uCacheCount);

	m_Ids = myalloc(uint, m_uCacheCount);

	m_dEdgeLength1 = myalloc(double, m_uCacheCount);
	m_dEdgeLength2 = myalloc(double, m_uCacheCount);
	m_dEdgeLength3 = myalloc(double, m_uCacheCount);
	m_dHeight = myalloc(double, m_uCacheCount);

	m_bHasEdgeLength1 = myalloc(bool, m_uCacheCount);
	m_bHasEdgeLength2 = myalloc(bool, m_uCacheCount);
	m_bHasEdgeLength3 = myalloc(bool, m_uCacheCount);
	m_bHasHeight = myalloc(bool, m_uCacheCount);

	m_ptrName = myalloc(char *, m_uCacheCount);

	for (uint uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		m_uNeighbor1[uNodeIndex] = NULL_NEIGHBOR;
		m_uNeighbor2[uNodeIndex] = NULL_NEIGHBOR;
		m_uNeighbor3[uNodeIndex] = NULL_NEIGHBOR;
		m_bHasEdgeLength1[uNodeIndex] = false;
		m_bHasEdgeLength2[uNodeIndex] = false;
		m_bHasEdgeLength3[uNodeIndex] = false;
		m_bHasHeight[uNodeIndex] = false;
		m_dEdgeLength1[uNodeIndex] = DBL_MAX;
		m_dEdgeLength2[uNodeIndex] = DBL_MAX;
		m_dEdgeLength3[uNodeIndex] = DBL_MAX;
		m_dHeight[uNodeIndex] = DBL_MAX;
		m_ptrName[uNodeIndex] = 0;
		m_Ids[uNodeIndex] = UINT_MAX;
		}
	}

void Tree::AssertAreNeighbors(uint uNodeIndex1, uint uNodeIndex2) const
	{
	if (uNodeIndex1 >= m_uNodeCount || uNodeIndex2 >= m_uNodeCount)
		Die("AssertAreNeighbors(%u,%u), are %u nodes",
		  uNodeIndex1, uNodeIndex2, m_uNodeCount);

	if (m_uNeighbor1[uNodeIndex1] != uNodeIndex2 &&
	  m_uNeighbor2[uNodeIndex1] != uNodeIndex2 &&
	  m_uNeighbor3[uNodeIndex1] != uNodeIndex2)
		{
		LogMe();
		Die("AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
		}

	if (m_uNeighbor1[uNodeIndex2] != uNodeIndex1 &&
	  m_uNeighbor2[uNodeIndex2] != uNodeIndex1 &&
	  m_uNeighbor3[uNodeIndex2] != uNodeIndex1)
		{
		LogMe();
		Die("AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
		}

	bool Has12 = HasEdgeLength(uNodeIndex1, uNodeIndex2);
	bool Has21 = HasEdgeLength(uNodeIndex2, uNodeIndex1);
	if (Has12 != Has21)
		{
		HasEdgeLength(uNodeIndex1, uNodeIndex2);
		HasEdgeLength(uNodeIndex2, uNodeIndex1);
		LogMe();
		Log("HasEdgeLength(%u, %u)=%c HasEdgeLength(%u, %u)=%c\n",
		  uNodeIndex1,
		  uNodeIndex2,
		  Has12 ? 'T' : 'F',
		  uNodeIndex2,
		  uNodeIndex1,
		  Has21 ? 'T' : 'F');

		Die("Tree::AssertAreNeighbors, HasEdgeLength not symmetric");
		}

	if (Has12)
		{
		double d12 = GetEdgeLength(uNodeIndex1, uNodeIndex2);
		double d21 = GetEdgeLength(uNodeIndex2, uNodeIndex1);
		if (d12 != d21)
			{
			LogMe();
			Die("Tree::AssertAreNeighbors, Edge length disagrees %u-%u=%.3g, %u-%u=%.3g",
			  uNodeIndex1, uNodeIndex2, d12,
			  uNodeIndex2, uNodeIndex1, d21);
			}
		}
	}

void Tree::ValidateNode(uint uNodeIndex) const
	{
	if (uNodeIndex >= m_uNodeCount)
		Die("ValidateNode(%u), %u nodes", uNodeIndex, m_uNodeCount);

	const uint uNeighborCount = GetNeighborCount(uNodeIndex);

	if (2 == uNeighborCount)
		{
		if (!m_bRooted)
			{
			LogMe();
			Die("Tree::ValidateNode: Node %u has two neighbors, tree is not rooted",
			 uNodeIndex);
			}
		if (uNodeIndex != m_uRootNodeIndex)
			{
			LogMe();
			Die("Tree::ValidateNode: Node %u has two neighbors, but not root node=%u",
			 uNodeIndex, m_uRootNodeIndex);
			}
		}

	const uint n1 = m_uNeighbor1[uNodeIndex];
	const uint n2 = m_uNeighbor2[uNodeIndex];
	const uint n3 = m_uNeighbor3[uNodeIndex];

	if (NULL_NEIGHBOR == n2 && NULL_NEIGHBOR != n3)
		{
		LogMe();
		Die("Tree::ValidateNode, n2=null, n3!=null", uNodeIndex);
		}
	if (NULL_NEIGHBOR == n3 && NULL_NEIGHBOR != n2)
		{
		LogMe();
		Die("Tree::ValidateNode, n3=null, n2!=null", uNodeIndex);
		}

	if (n1 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n1);
	if (n2 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n2);
	if (n3 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n3);

	if (n1 != NULL_NEIGHBOR && (n1 == n2 || n1 == n3))
		{
		LogMe();
		Die("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}
	if (n2 != NULL_NEIGHBOR && (n2 == n1 || n2 == n3))
		{
		LogMe();
		Die("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}
	if (n3 != NULL_NEIGHBOR && (n3 == n1 || n3 == n2))
		{
		LogMe();
		Die("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}

	if (IsRooted())
		{
		if (NULL_NEIGHBOR == GetParent(uNodeIndex))
			{
			if (uNodeIndex != m_uRootNodeIndex)
				{
				LogMe();
				Die("Tree::ValiateNode(%u), no parent", uNodeIndex);
				}
			}
		else if (GetLeft(GetParent(uNodeIndex)) != uNodeIndex &&
		  GetRight(GetParent(uNodeIndex)) != uNodeIndex)
			{
			LogMe();
			Die("Tree::ValidateNode(%u), parent / child mismatch", uNodeIndex);
			}
		}
	}

void Tree::Validate() const
	{
	for (uint uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		ValidateNode(uNodeIndex);
	}

uint Tree::GetSibling(uint Node) const
	{
	if (Node == UINT_MAX)
		return UINT_MAX;
	uint Parent = GetParent(Node);
	if (Parent == UINT_MAX)
		return UINT_MAX;
	uint ParentLeft = GetLeft(Parent);
	uint ParentRight = GetRight(Parent);
	asserta(ParentLeft != UINT_MAX);
	asserta(ParentRight != UINT_MAX);
	if (ParentLeft == Node)
		return ParentRight;
	if (ParentRight == Node)
		return ParentLeft;
	asserta(false);
	return UINT_MAX;
	}

bool Tree::IsEdge(uint uNodeIndex1, uint uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);

	return m_uNeighbor1[uNodeIndex1] == uNodeIndex2 ||
	  m_uNeighbor2[uNodeIndex1] == uNodeIndex2 ||
	  m_uNeighbor3[uNodeIndex1] == uNodeIndex2;
	}

double Tree::GetEdgeLength(uint uNodeIndex1, uint uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);
	if (!HasEdgeLength(uNodeIndex1, uNodeIndex2))
		{
		LogMe();
		Die("Missing edge length in tree %u-%u", uNodeIndex1, uNodeIndex2);
		}

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		return m_dEdgeLength1[uNodeIndex1];
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		return m_dEdgeLength2[uNodeIndex1];
	assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
	return m_dEdgeLength3[uNodeIndex1];
	}

void Tree::ExpandCache()
	{
	const uint uNodeCount = 100;
	uint uNewCacheCount = m_uCacheCount + uNodeCount;
	uint *uNewNeighbor1 = myalloc(uint, uNewCacheCount);
	uint *uNewNeighbor2 = myalloc(uint, uNewCacheCount);
	uint *uNewNeighbor3 = myalloc(uint, uNewCacheCount);

	uint *uNewIds = myalloc(uint, uNewCacheCount);
	memset(uNewIds, 0xff, uNewCacheCount*sizeof(uint));

	double *dNewEdgeLength1 = myalloc(double, uNewCacheCount);
	double *dNewEdgeLength2 = myalloc(double, uNewCacheCount);
	double *dNewEdgeLength3 = myalloc(double, uNewCacheCount);
	double *dNewHeight = myalloc(double, uNewCacheCount);

	bool *bNewHasEdgeLength1 = myalloc(bool, uNewCacheCount);
	bool *bNewHasEdgeLength2 = myalloc(bool, uNewCacheCount);
	bool *bNewHasEdgeLength3 = myalloc(bool, uNewCacheCount);
	bool *bNewHasHeight = myalloc(bool, uNewCacheCount);

	char **ptrNewName = myalloc(char *, uNewCacheCount);
	memset(ptrNewName, 0, uNewCacheCount*sizeof(char *));

	if (m_uCacheCount > 0)
		{
		const uint uUnsignedBytes = m_uCacheCount*sizeof(uint);
		memcpy(uNewNeighbor1, m_uNeighbor1, uUnsignedBytes);
		memcpy(uNewNeighbor2, m_uNeighbor2, uUnsignedBytes);
		memcpy(uNewNeighbor3, m_uNeighbor3, uUnsignedBytes);

		memcpy(uNewIds, m_Ids, uUnsignedBytes);

		const uint uEdgeBytes = m_uCacheCount*sizeof(double);
		memcpy(dNewEdgeLength1, m_dEdgeLength1, uEdgeBytes);
		memcpy(dNewEdgeLength2, m_dEdgeLength2, uEdgeBytes);
		memcpy(dNewEdgeLength3, m_dEdgeLength3, uEdgeBytes);
		memcpy(dNewHeight, m_dHeight, uEdgeBytes);

		const uint uBoolBytes = m_uCacheCount*sizeof(bool);
		memcpy(bNewHasEdgeLength1, m_bHasEdgeLength1, uBoolBytes);
		memcpy(bNewHasEdgeLength2, m_bHasEdgeLength2, uBoolBytes);
		memcpy(bNewHasEdgeLength3, m_bHasEdgeLength3, uBoolBytes);
		memcpy(bNewHasHeight, m_bHasHeight, uBoolBytes);

		const uint uNameBytes = m_uCacheCount*sizeof(char *);
		memcpy(ptrNewName, m_ptrName, uNameBytes);

		myfree(m_uNeighbor1);
		myfree(m_uNeighbor2);
		myfree(m_uNeighbor3);

		myfree(m_Ids);

		myfree(m_dEdgeLength1);
		myfree(m_dEdgeLength2);
		myfree(m_dEdgeLength3);

		myfree(m_bHasEdgeLength1);
		myfree(m_bHasEdgeLength2);
		myfree(m_bHasEdgeLength3);
		myfree(m_bHasHeight);

		myfree(m_ptrName);
		}
	m_uCacheCount = uNewCacheCount;
	m_uNeighbor1 = uNewNeighbor1;
	m_uNeighbor2 = uNewNeighbor2;
	m_uNeighbor3 = uNewNeighbor3;
	m_Ids = uNewIds;
	m_dEdgeLength1 = dNewEdgeLength1;
	m_dEdgeLength2 = dNewEdgeLength2;
	m_dEdgeLength3 = dNewEdgeLength3;
	m_dHeight = dNewHeight;
	m_bHasEdgeLength1 = bNewHasEdgeLength1;
	m_bHasEdgeLength2 = bNewHasEdgeLength2;
	m_bHasEdgeLength3 = bNewHasEdgeLength3;
	m_bHasHeight = bNewHasHeight;
	m_ptrName = ptrNewName;
	}

// Creates tree with single node, no edges.
// Root node always has index 0.
void Tree::CreateRooted()
	{
	Clear();
	ExpandCache();
	m_uNodeCount = 1;

	m_uNeighbor1[0] = NULL_NEIGHBOR;
	m_uNeighbor2[0] = NULL_NEIGHBOR;
	m_uNeighbor3[0] = NULL_NEIGHBOR;

	m_bHasEdgeLength1[0] = false;
	m_bHasEdgeLength2[0] = false;
	m_bHasEdgeLength3[0] = false;
	m_bHasHeight[0] = false;

	m_uRootNodeIndex = 0;
	m_bRooted = true;

#if	DEBUG
	Validate();
#endif
	}

// Creates unrooted tree with single edge.
// Nodes for that edge are always 0 and 1.
void Tree::CreateUnrooted(double dEdgeLength)
	{
	Clear();
	ExpandCache();

	m_uNeighbor1[0] = 1;
	m_uNeighbor2[0] = NULL_NEIGHBOR;
	m_uNeighbor3[0] = NULL_NEIGHBOR;

	m_uNeighbor1[1] = 0;
	m_uNeighbor2[1] = NULL_NEIGHBOR;
	m_uNeighbor3[1] = NULL_NEIGHBOR;

	m_dEdgeLength1[0] = dEdgeLength;
	m_dEdgeLength1[1] = dEdgeLength;

	m_bHasEdgeLength1[0] = true;
	m_bHasEdgeLength1[1] = true;

	m_bRooted = false;

#if	DEBUG
	Validate();
#endif
	}

void Tree::SetLeafName(uint uNodeIndex, const char *ptrName)
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	free(m_ptrName[uNodeIndex]);
	m_ptrName[uNodeIndex] = mystrsave(ptrName);
	}

void Tree::SetLeafId(uint uNodeIndex, uint uId)
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	m_Ids[uNodeIndex] = uId;
	}

const char *Tree::GetLeafName(uint uNodeIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	return m_ptrName[uNodeIndex];
	}

uint Tree::GetLeafId(uint uNodeIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	return m_Ids[uNodeIndex];
	}

// Append a new branch.
// This adds two new nodes and joins them to an existing leaf node.
// Return value is k, new nodes have indexes k and k+1 respectively.
uint Tree::AppendBranch(uint uExistingLeafIndex)
	{
	if (0 == m_uNodeCount)
		Die("Tree::AppendBranch: tree has not been created");

#if	DEBUG
	assert(uExistingLeafIndex < m_uNodeCount);
	if (!IsLeaf(uExistingLeafIndex))
		{
		LogMe();
		Die("AppendBranch(%u): not leaf", uExistingLeafIndex);
		}
#endif

	if (m_uNodeCount >= m_uCacheCount - 2)
		ExpandCache();

	const uint uNewLeaf1 = m_uNodeCount;
	const uint uNewLeaf2 = m_uNodeCount + 1;

	m_uNodeCount += 2;

	assert(m_uNeighbor2[uExistingLeafIndex] == NULL_NEIGHBOR);
	assert(m_uNeighbor3[uExistingLeafIndex] == NULL_NEIGHBOR);

	m_uNeighbor2[uExistingLeafIndex] = uNewLeaf1;
	m_uNeighbor3[uExistingLeafIndex] = uNewLeaf2;

	m_uNeighbor1[uNewLeaf1] = uExistingLeafIndex;
	m_uNeighbor1[uNewLeaf2] = uExistingLeafIndex;

	m_uNeighbor2[uNewLeaf1] = NULL_NEIGHBOR;
	m_uNeighbor2[uNewLeaf2] = NULL_NEIGHBOR;

	m_uNeighbor3[uNewLeaf1] = NULL_NEIGHBOR;
	m_uNeighbor3[uNewLeaf2] = NULL_NEIGHBOR;

	m_dEdgeLength2[uExistingLeafIndex] = 0;
	m_dEdgeLength3[uExistingLeafIndex] = 0;

	m_dEdgeLength1[uNewLeaf1] = 0;
	m_dEdgeLength2[uNewLeaf1] = 0;
	m_dEdgeLength3[uNewLeaf1] = 0;

	m_dEdgeLength1[uNewLeaf2] = 0;
	m_dEdgeLength2[uNewLeaf2] = 0;
	m_dEdgeLength3[uNewLeaf2] = 0;

	m_bHasEdgeLength1[uNewLeaf1] = false;
	m_bHasEdgeLength2[uNewLeaf1] = false;
	m_bHasEdgeLength3[uNewLeaf1] = false;

	m_bHasEdgeLength1[uNewLeaf2] = false;
	m_bHasEdgeLength2[uNewLeaf2] = false;
	m_bHasEdgeLength3[uNewLeaf2] = false;

	m_bHasHeight[uNewLeaf1] = false;
	m_bHasHeight[uNewLeaf2] = false;

	m_Ids[uNewLeaf1] = UINT_MAX;
	m_Ids[uNewLeaf2] = UINT_MAX;
	return uNewLeaf1;
	}

void Tree::LogMe() const
	{
	Log("Tree::LogMe %u nodes, ", m_uNodeCount);

	if (IsRooted())
		{
		Log("rooted.\n");
		Log("\n");
		Log("Index  Parnt  LengthP  Left   LengthL  Right  LengthR     Id  Name\n");
		Log("-----  -----  -------  ----   -------  -----  -------  -----  ----\n");
		}
	else
		{
		Log("unrooted.\n");
		Log("\n");
		Log("Index  Nbr_1  Length1  Nbr_2  Length2  Nbr_3  Length3     Id  Name\n");
		Log("-----  -----  -------  -----  -------  -----  -------  -----  ----\n");
		}

	for (uint uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		Log("%5u  ", uNodeIndex);
		const uint n1 = m_uNeighbor1[uNodeIndex];
		const uint n2 = m_uNeighbor2[uNodeIndex];
		const uint n3 = m_uNeighbor3[uNodeIndex];
		if (NULL_NEIGHBOR != n1)
			{
			Log("%5u  ", n1);
			if (m_bHasEdgeLength1[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength1[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (NULL_NEIGHBOR != n2)
			{
			Log("%5u  ", n2);
			if (m_bHasEdgeLength2[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength2[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (NULL_NEIGHBOR != n3)
			{
			Log("%5u  ", n3);
			if (m_bHasEdgeLength3[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength3[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (m_Ids != 0 && IsLeaf(uNodeIndex))
			{
			uint uId = m_Ids[uNodeIndex];
			if (uId == UINT_MAX)
				Log("    *");
			else
				Log("%5u", uId);
			}
		else
			Log("     ");

		if (m_bRooted && uNodeIndex == m_uRootNodeIndex)
			Log("  [ROOT] ");
		const char *ptrName = m_ptrName[uNodeIndex];
		if (ptrName != 0)
			Log("  %s", ptrName);
		Log("\n");
		}
	}

void Tree::SetEdgeLength(uint uNodeIndex1, uint uNodeIndex2,
  double dLength)
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);
	assert(IsEdge(uNodeIndex1, uNodeIndex2));

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		{
		m_dEdgeLength1[uNodeIndex1] = dLength;
		m_bHasEdgeLength1[uNodeIndex1] = true;
		}
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		{
		m_dEdgeLength2[uNodeIndex1] = dLength;
		m_bHasEdgeLength2[uNodeIndex1] = true;

		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
		m_dEdgeLength3[uNodeIndex1] = dLength;
		m_bHasEdgeLength3[uNodeIndex1] = true;
		}

	if (m_uNeighbor1[uNodeIndex2] == uNodeIndex1)
		{
		m_dEdgeLength1[uNodeIndex2] = dLength;
		m_bHasEdgeLength1[uNodeIndex2] = true;
		}
	else if (m_uNeighbor2[uNodeIndex2] == uNodeIndex1)
		{
		m_dEdgeLength2[uNodeIndex2] = dLength;
		m_bHasEdgeLength2[uNodeIndex2] = true;
		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex2] == uNodeIndex1);
		m_dEdgeLength3[uNodeIndex2] = dLength;
		m_bHasEdgeLength3[uNodeIndex2] = true;
		}
	}

uint Tree::UnrootFromFile()
	{
#if	TRACE
	Log("Before unroot:\n");
	LogMe();
#endif

	if (!m_bRooted)
		Die("Tree::Unroot, not rooted");

// Convention: root node is always node zero
	assert(IsRoot(0));
	assert(NULL_NEIGHBOR == m_uNeighbor1[0]);

	const uint uThirdNode = m_uNodeCount++;

	m_uNeighbor1[0] = uThirdNode;
	m_uNeighbor1[uThirdNode] = 0;

	m_uNeighbor2[uThirdNode] = NULL_NEIGHBOR;
	m_uNeighbor3[uThirdNode] = NULL_NEIGHBOR;

	m_dEdgeLength1[0] = 0;
	m_dEdgeLength1[uThirdNode] = 0;
	m_bHasEdgeLength1[uThirdNode] = true;

	m_bRooted = false;

#if	TRACE
	Log("After unroot:\n");
	LogMe();
#endif

	return uThirdNode;
	}

// In an unrooted tree, equivalent of GetLeft/Right is 
// GetFirst/SecondNeighbor.
// uNeighborIndex must be a known neighbor of uNodeIndex.
// This is the way to find the other two neighbor nodes of
// an internal node.
// The labeling as "First" and "Second" neighbor is arbitrary.
// Calling these functions on a leaf returns NULL_NEIGHBOR, as
// for GetLeft/Right.
uint Tree::GetFirstNeighbor(uint uNodeIndex, uint uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	assert(IsEdge(uNodeIndex, uNeighborIndex));

	for (uint n = 0; n < 3; ++n)
		{
		uint uNeighbor = GetNeighbor(uNodeIndex, n);
		if (NULL_NEIGHBOR != uNeighbor && uNeighborIndex != uNeighbor)
			return uNeighbor;
		}
	return NULL_NEIGHBOR;
	}

uint Tree::GetSecondNeighbor(uint uNodeIndex, uint uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	assert(IsEdge(uNodeIndex, uNeighborIndex));

	bool bFoundOne = false;
	for (uint n = 0; n < 3; ++n)
		{
		uint uNeighbor = GetNeighbor(uNodeIndex, n);
		if (NULL_NEIGHBOR != uNeighbor && uNeighborIndex != uNeighbor)
			{
			if (bFoundOne)
				return uNeighbor;
			else
				bFoundOne = true;
			}
		}
	return NULL_NEIGHBOR;
	}

// Compute the number of leaves in the sub-tree defined by an edge
// in an unrooted tree. Conceptually, the tree is cut at this edge,
// and uNodeIndex2 considered the root of the sub-tree.
uint Tree::GetLeafCountUnrooted(uint uNodeIndex1, uint uNodeIndex2,
  double *ptrdTotalDistance) const
	{
	assert(!IsRooted());

	if (IsLeaf(uNodeIndex2))
		{
		*ptrdTotalDistance = GetEdgeLength(uNodeIndex1, uNodeIndex2);
		return 1;
		}

// Recurse down the rooted sub-tree defined by cutting the edge
// and considering uNodeIndex2 as the root.
	const uint uLeft = GetFirstNeighbor(uNodeIndex2, uNodeIndex1);
	const uint uRight = GetSecondNeighbor(uNodeIndex2, uNodeIndex1);

	double dLeftDistance;
	double dRightDistance;

	const uint uLeftCount = GetLeafCountUnrooted(uNodeIndex2, uLeft,
	  &dLeftDistance);
	const uint uRightCount = GetLeafCountUnrooted(uNodeIndex2, uRight,
	  &dRightDistance);

	*ptrdTotalDistance = dLeftDistance + dRightDistance;
	return uLeftCount + uRightCount;
	}

bool Tree::HasEdgeLength(uint uNodeIndex1, uint uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount);
	assert(uNodeIndex2 < m_uNodeCount);
	assert(IsEdge(uNodeIndex1, uNodeIndex2));

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		return m_bHasEdgeLength1[uNodeIndex1];
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		return m_bHasEdgeLength2[uNodeIndex1];
	assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
	return m_bHasEdgeLength3[uNodeIndex1];
	}

void Tree::OrientParent(uint uNodeIndex, uint uParentNodeIndex)
	{
	if (NULL_NEIGHBOR == uNodeIndex)
		return;

	if (m_uNeighbor1[uNodeIndex] == uParentNodeIndex)
		;
	else if (m_uNeighbor2[uNodeIndex] == uParentNodeIndex)
		{
		double dEdgeLength2 = m_dEdgeLength2[uNodeIndex];
		m_uNeighbor2[uNodeIndex] = m_uNeighbor1[uNodeIndex];
		m_dEdgeLength2[uNodeIndex] = m_dEdgeLength1[uNodeIndex];
		m_uNeighbor1[uNodeIndex] = uParentNodeIndex;
		m_dEdgeLength1[uNodeIndex] = dEdgeLength2;
		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex] == uParentNodeIndex);
		double dEdgeLength3 = m_dEdgeLength3[uNodeIndex];
		m_uNeighbor3[uNodeIndex] = m_uNeighbor1[uNodeIndex];
		m_dEdgeLength3[uNodeIndex] = m_dEdgeLength1[uNodeIndex];
		m_uNeighbor1[uNodeIndex] = uParentNodeIndex;
		m_dEdgeLength1[uNodeIndex] = dEdgeLength3;
		}

	OrientParent(m_uNeighbor2[uNodeIndex], uNodeIndex);
	OrientParent(m_uNeighbor3[uNodeIndex], uNodeIndex);
	}

uint Tree::FirstDepthFirstNode() const
	{
	assert(IsRooted());

// Descend via left branches until we hit a leaf
	uint uNodeIndex = m_uRootNodeIndex;
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetLeft(uNodeIndex);
	return uNodeIndex;
	}

uint Tree::FirstDepthFirstNodeR() const
	{
	assert(IsRooted());

// Descend via left branches until we hit a leaf
	uint uNodeIndex = m_uRootNodeIndex;
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetRight(uNodeIndex);
	return uNodeIndex;
	}

uint Tree::NextDepthFirstNode(uint uNodeIndex) const
	{
#if	TRACE
	Log("NextDepthFirstNode(%3u) ", uNodeIndex);
#endif

	assert(IsRooted());
	assert(uNodeIndex < m_uNodeCount);

	if (IsRoot(uNodeIndex))
		{
#if	TRACE
		Log(">> Node %u is root, end of traversal\n", uNodeIndex);
#endif
		return NULL_NEIGHBOR;
		}

	uint uParent = GetParent(uNodeIndex);
	if (GetRight(uParent) == uNodeIndex)
		{
#if	TRACE
		Log(">> Is right branch, return parent=%u\n", uParent);
#endif
		return uParent;
		}

	uNodeIndex = GetRight(uParent);
#if	TRACE
		Log(">> Descend left from right sibling=%u ... ", uNodeIndex);
#endif
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetLeft(uNodeIndex);

#if	TRACE
	Log("bottom out at leaf=%u\n", uNodeIndex);
#endif
	return uNodeIndex;
	}

uint Tree::NextDepthFirstNodeR(uint uNodeIndex) const
	{
#if	TRACE
	Log("NextDepthFirstNode(%3u) ", uNodeIndex);
#endif

	assert(IsRooted());
	assert(uNodeIndex < m_uNodeCount);

	if (IsRoot(uNodeIndex))
		{
#if	TRACE
		Log(">> Node %u is root, end of traversal\n", uNodeIndex);
#endif
		return NULL_NEIGHBOR;
		}

	uint uParent = GetParent(uNodeIndex);
	if (GetLeft(uParent) == uNodeIndex)
		{
#if	TRACE
		Log(">> Is left branch, return parent=%u\n", uParent);
#endif
		return uParent;
		}

	uNodeIndex = GetLeft(uParent);
#if	TRACE
		Log(">> Descend right from left sibling=%u ... ", uNodeIndex);
#endif
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetRight(uNodeIndex);

#if	TRACE
	Log("bottom out at leaf=%u\n", uNodeIndex);
#endif
	return uNodeIndex;
	}

static void GetMaxString(const vector<string> &v, string &MaxStr)
	{
	asserta(!v.empty());
	MaxStr = v[0];
	for (uint i = 1; i < SIZE(v); ++i)
		MaxStr = max(MaxStr, v[i]);
	}

static bool CompareLabels(const vector<string> &Labels1,
  const vector<string> &Labels2)
	{
	string Max1;
	string Max2;
	GetMaxString(Labels1, Max1);
	GetMaxString(Labels2, Max2);
	bool Gt = (Max1 > Max2);
	return Gt;
	}

uint Tree::Ladderize(bool MoreRight)
	{
	const uint NodeCount = GetNodeCount();
	uint RotatedCount = 0;
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (IsLeaf(Node))
			continue;

		uint Left = GetLeft(Node);
		uint Right = GetRight(Node);
		uint NLeft = GetSubtreeLeafCount(Left);
		uint NRight = GetSubtreeLeafCount(Right);

		bool DoRotate = (MoreRight ? NRight < NLeft : NLeft < NRight);
		if (NLeft == NRight)
			{
			vector<string> LeftLabels;
			vector<string> RightLabels;
			GetSubtreeLeafLabels(Left, LeftLabels);
			GetSubtreeLeafLabels(Right, RightLabels);
			DoRotate = CompareLabels(LeftLabels, RightLabels);
			}
		if (DoRotate)
			{
			++RotatedCount;
			uint Left = GetLeft(Node);
			uint Right = GetRight(Node);
			m_uNeighbor2[Node] = Right;
			m_uNeighbor3[Node] = Left;
			}
		}
	return RotatedCount;
	}

void Tree::UnrootByDeletingRoot()
	{
	assert(IsRooted());
	assert(m_uNodeCount >= 3);

	const uint uLeft = GetLeft(m_uRootNodeIndex);
	const uint uRight = GetRight(m_uRootNodeIndex);

	m_uNeighbor1[uLeft] = uRight;
	m_uNeighbor1[uRight] = uLeft;

	bool bHasEdgeLength = HasEdgeLength(m_uRootNodeIndex, uLeft) &&
	  HasEdgeLength(m_uRootNodeIndex, uRight);
	if (bHasEdgeLength)
		{
		double dEdgeLength = GetEdgeLength(m_uRootNodeIndex, uLeft) +
		  GetEdgeLength(m_uRootNodeIndex, uRight);
		m_dEdgeLength1[uLeft] = dEdgeLength;
		m_dEdgeLength1[uRight] = dEdgeLength;
		}

// Remove root node entry from arrays
	const uint uMoveCount = m_uNodeCount - m_uRootNodeIndex;
	const uint uUnsBytes = uMoveCount*sizeof(uint);
	memmove(m_uNeighbor1 + m_uRootNodeIndex, m_uNeighbor1 + m_uRootNodeIndex + 1,
	  uUnsBytes);
	memmove(m_uNeighbor2 + m_uRootNodeIndex, m_uNeighbor2 + m_uRootNodeIndex + 1,
	  uUnsBytes);
	memmove(m_uNeighbor3 + m_uRootNodeIndex, m_uNeighbor3 + m_uRootNodeIndex + 1,
	  uUnsBytes);

	const uint uDoubleBytes = uMoveCount*sizeof(double);
	memmove(m_dEdgeLength1 + m_uRootNodeIndex, m_dEdgeLength1 + m_uRootNodeIndex + 1,
	  uDoubleBytes);
	memmove(m_dEdgeLength2 + m_uRootNodeIndex, m_dEdgeLength2 + m_uRootNodeIndex + 1,
	  uDoubleBytes);
	memmove(m_dEdgeLength3 + m_uRootNodeIndex, m_dEdgeLength3 + m_uRootNodeIndex + 1,
	  uDoubleBytes);

	const uint uBoolBytes = uMoveCount*sizeof(bool);
	memmove(m_bHasEdgeLength1 + m_uRootNodeIndex, m_bHasEdgeLength1 + m_uRootNodeIndex + 1,
	  uBoolBytes);
	memmove(m_bHasEdgeLength2 + m_uRootNodeIndex, m_bHasEdgeLength2 + m_uRootNodeIndex + 1,
	  uBoolBytes);
	memmove(m_bHasEdgeLength3 + m_uRootNodeIndex, m_bHasEdgeLength3 + m_uRootNodeIndex + 1,
	  uBoolBytes);

	const uint uPtrBytes = uMoveCount*sizeof(char *);
	memmove(m_ptrName + m_uRootNodeIndex, m_ptrName + m_uRootNodeIndex + 1, uPtrBytes);

	--m_uNodeCount;
	m_bRooted = false;

// Fix up table entries
	for (uint uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
#define DEC(x)	if (x != NULL_NEIGHBOR && x > m_uRootNodeIndex) --x;
		DEC(m_uNeighbor1[uNodeIndex])
		DEC(m_uNeighbor2[uNodeIndex])
		DEC(m_uNeighbor3[uNodeIndex])
#undef	DEC
		}

	Validate();
	}

uint Tree::GetLeafParent(uint uNodeIndex) const
	{
	assert(IsLeaf(uNodeIndex));

	if (IsRooted())
		return GetParent(uNodeIndex);

	if (m_uNeighbor1[uNodeIndex] != NULL_NEIGHBOR)
		return m_uNeighbor1[uNodeIndex];
	if (m_uNeighbor2[uNodeIndex] != NULL_NEIGHBOR)
		return m_uNeighbor2[uNodeIndex];
	return m_uNeighbor3[uNodeIndex];
	}

uint Tree::GetLCA(uint Node1, uint Node2) const
	{
	vector<uint> Path1;
	vector<uint> Path2;
	GetPathToRoot(Node1, Path1);
	GetPathToRoot(Node2, Path2);
	const uint N1 = SIZE(Path1);
	const uint N2 = SIZE(Path2);
	for (uint i = 0; i < N1; ++i)
		{
		uint AncNode1 = Path1[i];
		for (uint j = 0; j < N2; ++j)
			{
			if (Path2[j] == AncNode1)
				return AncNode1;
			}
		}
	asserta(false);
	return UINT_MAX;
	}

void Tree::GetPathToRoot(uint Node, vector<uint> &Path) const
	{
	if (!IsRooted())
		Die("GetPathToRoot(), not rooted");
	const uint NodeCount = GetNodeCount();
	Path.resize(0);
	for (;;)
		{
		asserta(Node < NodeCount);
		Path.push_back(Node);
		asserta(SIZE(Path) <= NodeCount);
		if (IsRoot(Node))
			return;
		Node = GetParent(Node);
		}
	}

// AncNode must be on path from Node to root
double Tree::GetDistance(uint Node, uint AncNode) const
	{
	if (Node == m_uRootNodeIndex && AncNode == UINT_MAX)
		return 0;

	vector<uint> Path;
	GetPathToRoot(Node, Path);
	double Distance = 0;
	const uint n = SIZE(Path);
	asserta(Path[0] == Node);
	for (uint i = 0; i < n; ++i)
		{
		uint PathNode = Path[i];
		if (PathNode == AncNode)
			return Distance;

	// Node1 is parent
		if (m_bHasEdgeLength1[PathNode])
			{
			double Length = m_dEdgeLength1[PathNode];
			Distance += Length;
			}
		}
	Die("GetDistance, not ancestor");
	return 0;
	}

void Tree::GetSubtreeLeafLabels(uint Node, vector<string> &Labels) const
	{
	Labels.clear();
	vector<uint> Leaves;
	AppendLeaves(Node, Leaves);
	uint n = SIZE(Leaves);
	for (uint i = 0; i < n; ++i)
		{
		uint LeafNode = Leaves[i];
		const char *Name = m_ptrName[LeafNode];
		Labels.push_back(string(Name));
		}
	}

void Tree::GetSubtreeLeafNodes(uint Node, vector<uint> &LeafNodes) const
	{
	AppendLeaves(Node, LeafNodes);
	}

uint Tree::GetSubtreeLeafCount(uint Node) const
	{
	if (Node == UINT_MAX)
		return 0;
	vector<uint> Leaves;
	AppendLeaves(Node, Leaves);
	uint n = SIZE(Leaves);
	return n;
	}

void Tree::GetLeafLabels(vector<string> &Labels) const
	{
	Labels.clear();
	const uint NodeCount = GetNodeCount();
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (IsLeaf(Node))
			{
			string Label;
			GetLabel(Node, Label);
			Labels.push_back(Label);
			}
		}
	}

void Tree::AppendLeaves(uint Node, vector<uint> &Leaves) const
	{
	uint NodeCount = GetNodeCount();
	uint LeafCount = GetLeafCount();
	asserta(Node < NodeCount);
	asserta(SIZE(Leaves) < LeafCount);

	if (IsLeaf(Node))
		Leaves.push_back(Node);
	else
		{
		uint Edge2 = m_uNeighbor2[Node];
		uint Edge3 = m_uNeighbor3[Node];
		AppendLeaves(Edge2, Leaves);
		AppendLeaves(Edge3, Leaves);
		}
	}

// TODO: This is not efficient for large trees, should cache.
double Tree::GetNodeHeight(uint uNodeIndex) const
	{
	if (!IsRooted())
		Die("Tree::GetNodeHeight: undefined unless rooted tree");
	
	if (IsLeaf(uNodeIndex))
		return 0.0;

	if (m_bHasHeight[uNodeIndex])
		return m_dHeight[uNodeIndex];

	const uint uLeft = GetLeft(uNodeIndex);
	const uint uRight = GetRight(uNodeIndex);
	double dLeftLength = GetEdgeLength(uNodeIndex, uLeft);
	double dRightLength = GetEdgeLength(uNodeIndex, uRight);

	if (dLeftLength < 0)
		dLeftLength = 0;
	if (dRightLength < 0)
		dRightLength = 0;

	const double dLeftHeight = dLeftLength + GetNodeHeight(uLeft);
	const double dRightHeight = dRightLength + GetNodeHeight(uRight);
	const double dHeight = (dLeftHeight + dRightHeight)/2;
	m_bHasHeight[uNodeIndex] = true;
	m_dHeight[uNodeIndex] = dHeight;
	return dHeight;
	}

uint Tree::GetNeighborSubscript(uint uNodeIndex, uint uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	if (uNeighborIndex == m_uNeighbor1[uNodeIndex])
		return 0;
	if (uNeighborIndex == m_uNeighbor2[uNodeIndex])
		return 1;
	if (uNeighborIndex == m_uNeighbor3[uNodeIndex])
		return 2;
	return NULL_NEIGHBOR;
	}

uint Tree::GetNeighbor(uint uNodeIndex, uint uNeighborSubscript) const
	{
	switch (uNeighborSubscript)
		{
	case 0:
		return m_uNeighbor1[uNodeIndex];
	case 1:
		return m_uNeighbor2[uNodeIndex];
	case 2:
		return m_uNeighbor3[uNodeIndex];
		}
	Die("Tree::GetNeighbor, sub=%u", uNeighborSubscript);
	return NULL_NEIGHBOR;
	}

// TODO: check if this is a performance issue, could cache a lookup table
uint Tree::LeafIndexToNodeIndex(uint uLeafIndex) const
	{
	const uint uNodeCount = GetNodeCount();
	uint uLeafCount = 0;
	for (uint uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (IsLeaf(uNodeIndex))
			{
			if (uLeafCount == uLeafIndex)
				return uNodeIndex;
			else
				++uLeafCount;
			}
		}
	Die("LeafIndexToNodeIndex: out of range");
	return 0;
	}

uint Tree::GetNodeIndex(const string &Label) const
	{
	return GetNodeIndex(Label.c_str());
	}

uint Tree::GetNodeIndex(const char *ptrName) const
	{
	const uint uNodeCount = GetNodeCount();
	for (uint uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		const char *ptrLeafName = m_ptrName[uNodeIndex];
		if (ptrLeafName != 0 && 0 == strcmp(ptrName, ptrLeafName))
			return uNodeIndex;
		}
	Die("Tree::GetLeafNodeIndex, name not found");
	return 0;
	}

void Tree::Copy(const Tree &tree)
	{
	const uint uNodeCount = tree.GetNodeCount();
	InitCache(uNodeCount);

	m_uNodeCount = uNodeCount;

	const size_t UnsignedBytes = uNodeCount*sizeof(uint);
	const size_t DoubleBytes = uNodeCount*sizeof(double);
	const size_t BoolBytes = uNodeCount*sizeof(bool);

	memcpy(m_uNeighbor1, tree.m_uNeighbor1, UnsignedBytes);
	memcpy(m_uNeighbor2, tree.m_uNeighbor2, UnsignedBytes);
	memcpy(m_uNeighbor3, tree.m_uNeighbor3, UnsignedBytes);

	memcpy(m_Ids, tree.m_Ids, UnsignedBytes);

	memcpy(m_dEdgeLength1, tree.m_dEdgeLength1, DoubleBytes);
	memcpy(m_dEdgeLength2, tree.m_dEdgeLength2, DoubleBytes);
	memcpy(m_dEdgeLength3, tree.m_dEdgeLength3, DoubleBytes);

	memcpy(m_dHeight, tree.m_dHeight, DoubleBytes);

	memcpy(m_bHasEdgeLength1, tree.m_bHasEdgeLength1, BoolBytes);
	memcpy(m_bHasEdgeLength2, tree.m_bHasEdgeLength2, BoolBytes);
	memcpy(m_bHasEdgeLength3, tree.m_bHasEdgeLength3, BoolBytes);

	memcpy(m_bHasHeight, tree.m_bHasHeight, BoolBytes);

	m_uRootNodeIndex = tree.m_uRootNodeIndex;
	m_bRooted = tree.m_bRooted;

	for (uint uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		if (tree.IsLeaf(uNodeIndex))
			{
			const char *ptrName = tree.GetLeafName(uNodeIndex);
			m_ptrName[uNodeIndex] = mystrsave(ptrName);
			}
		else
			m_ptrName[uNodeIndex] = 0;
		}

#if	DEBUG
	Validate();
#endif
	}

void Tree::ToVectors(vector<string> &Labels, 
  vector<uint> &Parents, vector<float> &Lengths) const
	{
	asserta(IsRooted());

	Labels.clear();
	Parents.clear();
	Lengths.clear();

	const uint NodeCount = GetNodeCount();
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		string Label;
		GetLabel(Node, Label);
		uint Parent = GetParent(Node);
		float Length = 0;
		if (Parent != UINT_MAX)
			Length = (float) GetEdgeLength(Node, Parent);

		Labels.push_back(Label);
		Parents.push_back(Parent);
		Lengths.push_back(Length);
		}
	}

void Tree::FromVectors(const vector<string> &Labels, 
  const vector<uint> &Parents, const vector<float> &Lengths)
	{
	Clear();
	const uint NodeCount = SIZE(Labels);
	asserta(SIZE(Parents) == NodeCount);
	asserta(SIZE(Lengths) == NodeCount);

	vector<uint> Lefts(NodeCount, UINT_MAX);
	vector<uint> Rights(NodeCount, UINT_MAX);
	uint Root = UINT_MAX;
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		uint Parent = Parents[Node];
		if (Parent == UINT_MAX)
			{
			asserta(Root == UINT_MAX);
			Root = Node;
			continue;
			}
		asserta(Parent < NodeCount);
		if (Lefts[Parent] == UINT_MAX)
			Lefts[Parent] = Node;
		else if (Rights[Parent] == UINT_MAX)
			Rights[Parent] = Node;
		else
			Die("Tree::FromVectors(), invalid vector topology");
		}
	asserta(Root != UINT_MAX);

	vector<uint> LeafNodes;
	vector<uint> IntNodes;
	vector<uint> NodeToLeafIndex(NodeCount, UINT_MAX);
	vector<uint> NodeToIntIndex(NodeCount, UINT_MAX);
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (Lefts[Node] == UINT_MAX)
			{
			asserta(Rights[Node] == UINT_MAX);
			uint LeafIndex = SIZE(LeafNodes);
			NodeToLeafIndex[Node] = LeafIndex;
			LeafNodes.push_back(Node);
			}
		else
			{
			asserta(Rights[Node] != UINT_MAX);
			uint IntIndex = SIZE(IntNodes);
			NodeToIntIndex[Node] = IntIndex;
			IntNodes.push_back(Node);
			}
		}

	uint LeafCount = SIZE(LeafNodes);
	uint IntCount = SIZE(IntNodes);
	asserta(LeafCount == (NodeCount + 1)/2);
	asserta(IntCount = LeafCount - 1);

	m_uNodeCount = NodeCount;
	InitCache(NodeCount);

// Leaves
	for (uint i = 0; i < LeafCount; ++i)
		{
		uint Node = LeafNodes[i];
		asserta(Node < SIZE(Labels));
		const string &Label = Labels[Node];
		asserta(Label != "");
		m_Ids[i] = i;
		m_ptrName[i] = mystrsave(Label.c_str());
		}

// Internal ndoes
	for (uint i = 0; i < IntCount; ++i)
		{
		uint Node = IntNodes[i];
		asserta(Node < SIZE(Lefts));
		asserta(Node < SIZE(Rights));

		uint NewNode = LeafCount + i;

		uint Left = Lefts[Node];
		uint Right = Rights[Node];

		uint LeftIntIndex = NodeToIntIndex[Left];
		uint RightIntIndex = NodeToIntIndex[Right];

		uint LeftLeafIndex = NodeToLeafIndex[Left];
		uint RightLeafIndex = NodeToLeafIndex[Right];

		uint NewLeftNode = UINT_MAX;
		if (LeftIntIndex == UINT_MAX)
			{
			asserta(LeftLeafIndex < LeafCount);
			NewLeftNode = LeftLeafIndex;
			}
		else
			{
			asserta(LeftIntIndex < IntCount);
			NewLeftNode = LeafCount + LeftIntIndex;
			}

		uint NewRightNode = UINT_MAX;
		if (RightIntIndex == UINT_MAX)
			{
			asserta(RightLeafIndex < LeafCount);
			NewRightNode = RightLeafIndex;
			}
		else
			{
			asserta(RightIntIndex < IntCount);
			NewRightNode = LeafCount + RightIntIndex;
			}

		m_uNeighbor2[NewNode] = NewLeftNode;
		m_uNeighbor3[NewNode] = NewRightNode;

		m_uNeighbor1[NewLeftNode] = NewNode;
		m_uNeighbor1[NewRightNode] = NewNode;

		m_bHasEdgeLength1[NewLeftNode] = false;
		m_bHasEdgeLength1[NewRightNode] = false;

		float LeftLength = Lengths[Left];
		float RightLength = Lengths[Right];
		m_bHasEdgeLength2[NewNode] = false;
		m_bHasEdgeLength3[NewNode] = false;
		if (LeftLength != MISSING_LENGTH)
			{
			m_bHasEdgeLength2[NewNode] = true;
			m_dEdgeLength2[NewNode] = LeftLength;

			m_bHasEdgeLength1[NewLeftNode] = true;
			m_dEdgeLength1[NewLeftNode] = LeftLength;
			}

		if (RightLength != MISSING_LENGTH)
			{
			m_bHasEdgeLength3[NewNode] = true;
			m_dEdgeLength3[NewNode] = RightLength;

			m_bHasEdgeLength1[NewRightNode] = true;
			m_dEdgeLength1[NewRightNode] = RightLength;
			}
		}

	uint NewRootIntIndex = NodeToIntIndex[Root];
	asserta(NewRootIntIndex != UINT_MAX);

	m_bRooted = true;
	m_uRootNodeIndex = LeafCount + NewRootIntIndex;

	Validate();
	}

// Create rooted tree from a vector description.
// Node indexes are 0..N-1 for leaves, N..2N-2 for
// internal nodes.
// Vector subscripts are i-N and have values for
// internal nodes only, but those values are node
// indexes 0..2N-2. So e.g. if N=6 and Left[2]=1,
// this means that the third internal node (node index 8)
// has the second leaf (node index 1) as its left child.
// uRoot gives the vector subscript of the root, so add N
// to get the node index.
void Tree::Create(uint uLeafCount, uint uRoot, const uint Left[],
  const uint Right[], const float LeftLength[], const float RightLength[],
  const uint LeafIds[], char **LeafNames)
	{
	Clear();

	m_uNodeCount = 2*uLeafCount - 1;
	InitCache(m_uNodeCount);

	for (uint uNodeIndex = 0; uNodeIndex < uLeafCount; ++uNodeIndex)
		{
		m_Ids[uNodeIndex] = LeafIds[uNodeIndex];
		m_ptrName[uNodeIndex] = mystrsave(LeafNames[uNodeIndex]);
		}

	for (uint uNodeIndex = uLeafCount; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		uint v = uNodeIndex - uLeafCount;
		uint uLeft = Left[v];
		uint uRight = Right[v];
		float fLeft = LeftLength[v];
		float fRight = RightLength[v];

		m_uNeighbor2[uNodeIndex] = uLeft;
		m_uNeighbor3[uNodeIndex] = uRight;

		m_bHasEdgeLength2[uNodeIndex] = true;
		m_bHasEdgeLength3[uNodeIndex] = true;

		m_dEdgeLength2[uNodeIndex] = fLeft;
		m_dEdgeLength3[uNodeIndex] = fRight;

		m_uNeighbor1[uLeft] = uNodeIndex;
		m_uNeighbor1[uRight] = uNodeIndex;

		m_dEdgeLength1[uLeft] = fLeft;
		m_dEdgeLength1[uRight] = fRight;

		m_bHasEdgeLength1[uLeft] = true;
		m_bHasEdgeLength1[uRight] = true;
		}

	m_bRooted = true;
	m_uRootNodeIndex = uRoot + uLeafCount;

	Validate();
	}

void Tree::GetSubtreeSizes(vector<uint> &Sizes) const
	{
	asserta(IsRooted());
	const uint NodeCount = GetNodeCount();
	Sizes.clear();
	Sizes.resize(NodeCount, UINT_MAX);
	uint Node = FirstDepthFirstNode();
	do
		{
		asserta(Node < NodeCount);
		if (IsLeaf(Node))
			Sizes[Node] = 1;
		else
			{
			uint Left = GetLeft(Node);
			uint Right = GetRight(Node);
			asserta(Left < NodeCount);
			asserta(Right < NodeCount);
			uint SizeLeft = Sizes[Left];
			uint SizeRight = Sizes[Right];
			asserta(SizeLeft != UINT_MAX && SizeLeft > 0);
			asserta(SizeRight != UINT_MAX && SizeRight > 0);
			asserta(Sizes[Node] == UINT_MAX);
			Sizes[Node] = SizeLeft + SizeRight;
			}
		Node = NextDepthFirstNode(Node);
		}
	while (Node != UINT_MAX);
	}
