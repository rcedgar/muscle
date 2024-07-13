#ifndef tree_h
#define tree_h

#include <limits.h>

class Clust;

const uint NULL_NEIGHBOR = UINT_MAX;
const double MISSING_LENGTH = DBL_MAX;

enum NEWICK_TOKEN_TYPE
	{
	NTT_Unknown,

// Returned from Tree::GetToken:
	NTT_Lparen,
	NTT_Rparen,
	NTT_Colon,
	NTT_Comma,
	NTT_Semicolon,
	NTT_String,

// Following are never returned from Tree::GetToken:
	NTT_SingleQuotedString,
	NTT_DoubleQuotedString,
	NTT_Comment
	};

class Tree
	{
public:

// Yuck. Data is made public for the convenience of Tree::Copy.
// There has to be a better way.
public:
	uint m_uNodeCount;
	uint m_uCacheCount;
	uint *m_uNeighbor1;
	uint *m_uNeighbor2;
	uint *m_uNeighbor3;
	double *m_dEdgeLength1;
	double *m_dEdgeLength2;
	double *m_dEdgeLength3;
	double *m_dHeight;
	bool *m_bHasEdgeLength1;
	bool *m_bHasEdgeLength2;
	bool *m_bHasEdgeLength3;
	bool *m_bHasHeight;
	uint *m_Ids;
	char **m_ptrName;
	bool m_bRooted = false;
	uint m_uRootNodeIndex = UINT_MAX;

	Tree()
		{
		m_uNodeCount = 0;
		m_uCacheCount = 0;
		m_uNeighbor1 = 0;
		m_uNeighbor2 = 0;
		m_uNeighbor3 = 0;
		m_dEdgeLength1 = 0;
		m_dEdgeLength2 = 0;
		m_dEdgeLength3 = 0;
		m_dHeight = 0;
		m_bHasEdgeLength1 = 0;
		m_bHasEdgeLength2 = 0;
		m_bHasEdgeLength3 = 0;
		m_bHasHeight = 0;
		m_ptrName = 0;
		m_Ids = 0;
		m_bRooted = false;
		m_uRootNodeIndex = UINT_MAX;
		}

	virtual ~Tree()
		{
		Clear();
		}

	void Clear()
		{
		for (uint n = 0; n < m_uNodeCount; ++n)
			free(m_ptrName[n]);

		m_uNodeCount = 0;
		m_uCacheCount = 0;

#define del(x)	if (x != 0) delete[] x; x = 0;
		del(m_uNeighbor1)
		del(m_uNeighbor2)
		del(m_uNeighbor3)
		del(m_dEdgeLength1)
		del(m_dEdgeLength2)
		del(m_dEdgeLength3)
		del(m_bHasEdgeLength1)
		del(m_bHasEdgeLength2)
		del(m_bHasEdgeLength3)
		del(m_ptrName)
		del(m_Ids)
		del(m_bHasHeight)
		del(m_dHeight)
#undef del

		m_uNeighbor1 = 0;
		m_uNeighbor2 = 0;
		m_uNeighbor3 = 0;
		m_dEdgeLength1 = 0;
		m_dEdgeLength2 = 0;
		m_dEdgeLength3 = 0;
		m_ptrName = 0;
		m_Ids = 0;
		m_uRootNodeIndex = 0;
		m_bHasHeight = 0;
		m_dHeight = 0;

		m_bRooted = false;
		}

// Creation and manipulation
	void CreateRooted();
	void CreateUnrooted(double dEdgeLength);

	void FromFile(const string &FileName);
	void FromFile(TextFile &File);
	//void FromClust(Clust &C);

	void Copy(const Tree &tree);

	void Create(uint uLeafCount, uint uRoot, const uint Left[],
	  const uint Right[], const float LeftLength[], const float RightLength[],
	  const uint LeafIds[], char *LeafNames[]);
	void FromVectors(const vector<string> &Labels, 
	  const vector<uint> &Parents, const vector<float> &Lengths);
	void ToVectors(vector<string> &Labels, 
	  vector<uint> &Parents, vector<float> &Lengths) const;

	uint AppendBranch(uint uExistingNodeIndex);
	void SetLeafName(uint uNodeIndex, const char *ptrName);
	void SetLeafId(uint uNodeIndex, uint uId);
	void SetEdgeLength(uint uNodeIndex1, uint uNodeIndex2,
	  double dLength);

	//void RootUnrootedTree(uint uNodeIndex1, uint uNodeIndex2);
	void UnrootByDeletingRoot();
	uint Ladderize(bool Right);

// Saving to file
	void ToFile(TextFile &File) const;
	void ToFile(const string &FileName) const;

// Accessor functions
	uint GetNodeCount() const
		{
		return m_uNodeCount;
		}

	uint GetLeafCount() const
		{
		if (m_bRooted)
			{
			assert(m_uNodeCount%2 == 1);
			return (m_uNodeCount + 1)/2;
			}
		else
			{
			assert(m_uNodeCount%2 == 0);
			return (m_uNodeCount + 2)/2;
			}
		}

	uint GetRoot() const
		{
		return m_uRootNodeIndex;
		}

	uint GetNeighbor(uint uNodeIndex, uint uNeighborSubscript) const;

	uint GetNeighbor1(uint uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor1[uNodeIndex];
		}

	uint GetNeighbor2(uint uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor2[uNodeIndex];
		}

	uint GetNeighbor3(uint uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor3[uNodeIndex];
		}

	uint GetParent(uint uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor1[uNodeIndex];
		}

	bool IsRooted() const
		{
		return m_bRooted;
		}

	uint GetSibling(uint uNodeIndex) const;

	uint GetLeft(uint uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor2[uNodeIndex];
		}

	uint GetRight(uint uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor3[uNodeIndex];
		}

	const char *GetName(uint uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_ptrName[uNodeIndex];
		}

	const void GetLabel(uint uNodeIndex, string &Label) const
		{
		assert(uNodeIndex < m_uNodeCount);
		const char *Name = m_ptrName[uNodeIndex];
		if (Name == 0)
			Label = "";
		else
			Label = string(Name);
		}

	uint GetRootNodeIndex() const
		{
		assert(m_bRooted);
		return m_uRootNodeIndex;
		}

	uint GetNeighborCount(uint uNodeIndex) const
		{
		const uint n1 = m_uNeighbor1[uNodeIndex];
		const uint n2 = m_uNeighbor2[uNodeIndex];
		const uint n3 = m_uNeighbor3[uNodeIndex];
		return (NULL_NEIGHBOR != n1) + (NULL_NEIGHBOR != n2) + (NULL_NEIGHBOR != n3);
		}

	bool IsLeaf(uint uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		if (1 == m_uNodeCount)
			return true;
		return 1 == GetNeighborCount(uNodeIndex);
		}

	bool IsRoot(uint uNodeIndex) const
		{
		return IsRooted() && m_uRootNodeIndex == uNodeIndex;
		}

	uint GetLeafId(uint uNodeIndex) const;
	uint GetNodeIndex(const char *ptrName) const;
	uint GetNodeIndex(const string &Label) const;
	bool IsEdge(uint uNodeIndex1, uint uNodeIndex2) const;
	bool HasEdgeLength(uint uNodeIndex1, uint uNodeIndex2) const;
	double GetEdgeLength(uint uNodeIndex1, uint uNodeIndex2) const;
	const char *GetLeafName(uint uNodeIndex) const;
	uint GetNeighborSubscript(uint uNodeIndex, uint uNeighborIndex) const;
	double GetNodeHeight(uint uNodeIndex) const;

	void AppendLeaves(uint Node, vector<uint> &Leaves) const;
	uint GetSubtreeLeafCount(uint Node) const;
	void GetSubtreeLeafLabels(uint Node, vector<string> &Labels) const;
	void GetSubtreeLeafNodes(uint Node, vector<uint> &Nodes) const;
	void GetPathToRoot(uint Node, vector<uint> &Path) const;
	double GetDistance(uint Node, uint AncNode) const;
	void GetLeafLabels(vector<string> &Labels) const;

// Depth-first traversal
	uint FirstDepthFirstNode() const;
	uint NextDepthFirstNode(uint uNodeIndex) const;

	uint FirstDepthFirstNodeR() const;
	uint NextDepthFirstNodeR(uint uNodeIndex) const;

// Equivalent of GetLeft/Right in unrooted tree, works in rooted tree too.
	uint GetFirstNeighbor(uint uNodeIndex, uint uNeighborIndex) const;
	uint GetSecondNeighbor(uint uNodeIndex, uint uNeighborIndex) const;

// Getting parent node in unrooted tree defined iff leaf
	uint GetLeafParent(uint uNodeIndex) const;

// Misc
	const char *NTTStr(NEWICK_TOKEN_TYPE NTT) const;
	//void FindCenterByLongestSpan(uint *ptrNodeIndex1,
	//  uint *ptrNodeIndex2) const;
	void PruneTree(const Tree &tree, uint Subfams[],
	  uint uSubfamCount, const char *LabelPrefix,
	  vector<string> &Labels);
	uint LeafIndexToNodeIndex(uint uLeafIndex) const;

// Debugging & trouble-shooting support
	void Validate() const;
	void ValidateNode(uint uNodeIndex) const;
	void AssertAreNeighbors(uint uNodeIndex1, uint uNodeIndex2) const;
	void LogMe() const;
	uint GetLCA(uint Node1, uint Node2) const;
	void GetSubtreeSizes(vector<uint> &Sizes) const;

private:
	uint UnrootFromFile();
	NEWICK_TOKEN_TYPE GetTokenVerbose(TextFile &File, char szToken[],
	  uint uBytes) const
		{
		NEWICK_TOKEN_TYPE NTT = GetToken(File, szToken, uBytes);
		Log("GetToken %10.10s  %s\n", NTTStr(NTT), szToken);
		return NTT;
		}

	void InitCache(uint uCacheCount);
	void ExpandCache();
	NEWICK_TOKEN_TYPE GetToken(TextFile &File, char szToken[], uint uBytes) const;
	bool GetGroupFromFile(TextFile &File, uint uNodeIndex, double *ptrdEdgeLength);
	uint GetLeafCountUnrooted(uint uNodeIndex1, uint uNodeIndex2,
	  double *ptrdTotalDistance) const;
	void ToFileNodeRooted(TextFile &File, uint uNodeIndex) const;
	void ToFileNodeUnrooted(TextFile &File, uint uNodeIndex, uint uParent) const;
	void OrientParent(uint uNodeIndex, uint uParentNodeIndex);
	//double FromClustNode(const Clust &C, uint uClustNodeIndex, uint uPhyNodeIndex);
	uint GetAnyNonLeafNode() const;
	};

struct PhyEnumEdgeState
	{
	PhyEnumEdgeState()
		{
		m_bInit = false;
		m_uNodeIndex1 = NULL_NEIGHBOR;
		m_uNodeIndex2 = NULL_NEIGHBOR;
		}
	bool m_bInit;
	uint m_uNodeIndex1;
	uint m_uNodeIndex2;
	};

#endif // tree_h
