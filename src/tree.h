#ifndef tree_h
#define tree_h

#include <limits.h>

class Clust;

const unsigned NULL_NEIGHBOR = UINT_MAX;
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
		}
	virtual ~Tree()
		{
		Clear();
		}

	void Clear()
		{
		for (unsigned n = 0; n < m_uNodeCount; ++n)
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
	void FromClust(Clust &C);

	void Copy(const Tree &tree);

	void Create(unsigned uLeafCount, unsigned uRoot, const unsigned Left[],
	  const unsigned Right[], const float LeftLength[], const float RightLength[],
	  const unsigned LeafIds[], char *LeafNames[]);
	void FromVectors(const vector<string> &Labels, 
	  const vector<uint> &Parents, const vector<float> &Lengths);
	void ToVectors(vector<string> &Labels, 
	  vector<uint> &Parents, vector<float> &Lengths) const;

	unsigned AppendBranch(unsigned uExistingNodeIndex);
	void SetLeafName(unsigned uNodeIndex, const char *ptrName);
	void SetLeafId(unsigned uNodeIndex, unsigned uId);
	void SetEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2,
	  double dLength);

	void RootUnrootedTree(unsigned uNodeIndex1, unsigned uNodeIndex2);
	void UnrootByDeletingRoot();
	uint Ladderize(bool Right);

// Saving to file
	void ToFile(TextFile &File) const;
	void ToFile(const string &FileName) const;

// Accessor functions
	unsigned GetNodeCount() const
		{
		return m_uNodeCount;
		}

	unsigned GetLeafCount() const
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

	unsigned GetNeighbor(unsigned uNodeIndex, unsigned uNeighborSubscript) const;

	unsigned GetNeighbor1(unsigned uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor1[uNodeIndex];
		}

	unsigned GetNeighbor2(unsigned uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor2[uNodeIndex];
		}

	unsigned GetNeighbor3(unsigned uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		return m_uNeighbor3[uNodeIndex];
		}

	unsigned GetParent(unsigned uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor1[uNodeIndex];
		}

	bool IsRooted() const
		{
		return m_bRooted;
		}

	uint GetSibling(uint uNodeIndex) const;

	unsigned GetLeft(unsigned uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor2[uNodeIndex];
		}

	unsigned GetRight(unsigned uNodeIndex) const
		{
		assert(m_bRooted && uNodeIndex < m_uNodeCount);
		return m_uNeighbor3[uNodeIndex];
		}

	const char *GetName(unsigned uNodeIndex) const
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

	unsigned GetRootNodeIndex() const
		{
		assert(m_bRooted);
		return m_uRootNodeIndex;
		}

	unsigned GetNeighborCount(unsigned uNodeIndex) const
		{
		const unsigned n1 = m_uNeighbor1[uNodeIndex];
		const unsigned n2 = m_uNeighbor2[uNodeIndex];
		const unsigned n3 = m_uNeighbor3[uNodeIndex];
		return (NULL_NEIGHBOR != n1) + (NULL_NEIGHBOR != n2) + (NULL_NEIGHBOR != n3);
		}

	bool IsLeaf(unsigned uNodeIndex) const
		{
		assert(uNodeIndex < m_uNodeCount);
		if (1 == m_uNodeCount)
			return true;
		return 1 == GetNeighborCount(uNodeIndex);
		}

	bool IsRoot(unsigned uNodeIndex) const
		{
		return IsRooted() && m_uRootNodeIndex == uNodeIndex;
		}

	unsigned GetLeafId(unsigned uNodeIndex) const;
	unsigned GetNodeIndex(const char *ptrName) const;
	unsigned GetNodeIndex(const string &Label) const;
	bool IsEdge(unsigned uNodeIndex1, unsigned uNodeIndex2) const;
	bool HasEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2) const;
	double GetEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2) const;
	const char *GetLeafName(unsigned uNodeIndex) const;
	unsigned GetNeighborSubscript(unsigned uNodeIndex, unsigned uNeighborIndex) const;
	double GetNodeHeight(unsigned uNodeIndex) const;

	void AppendLeaves(uint Node, vector<uint> &Leaves) const;
	uint GetSubtreeLeafCount(uint Node) const;
	void GetSubtreeLeafLabels(uint Node, vector<string> &Labels) const;
	void GetSubtreeLeafNodes(uint Node, vector<uint> &Labels) const;
	void GetPathToRoot(uint Node, vector<uint> &Path) const;
	double GetDistance(uint Node, uint AncNode) const;
	void GetLeafLabels(vector<string> &Labels) const;

// Depth-first traversal
	unsigned FirstDepthFirstNode() const;
	unsigned NextDepthFirstNode(unsigned uNodeIndex) const;

	unsigned FirstDepthFirstNodeR() const;
	unsigned NextDepthFirstNodeR(unsigned uNodeIndex) const;

// Equivalent of GetLeft/Right in unrooted tree, works in rooted tree too.
	unsigned GetFirstNeighbor(unsigned uNodeIndex, unsigned uNeighborIndex) const;
	unsigned GetSecondNeighbor(unsigned uNodeIndex, unsigned uNeighborIndex) const;

// Getting parent node in unrooted tree defined iff leaf
	unsigned GetLeafParent(unsigned uNodeIndex) const;

// Misc
	const char *NTTStr(NEWICK_TOKEN_TYPE NTT) const;
	void FindCenterByLongestSpan(unsigned *ptrNodeIndex1,
	  unsigned *ptrNodeIndex2) const;
	void PruneTree(const Tree &tree, unsigned Subfams[],
	  unsigned uSubfamCount);
	unsigned LeafIndexToNodeIndex(unsigned uLeafIndex) const;

// Debugging & trouble-shooting support
	void Validate() const;
	void ValidateNode(unsigned uNodeIndex) const;
	void AssertAreNeighbors(unsigned uNodeIndex1, unsigned uNodeIndex2) const;
	void LogMe() const;
	uint GetLCA(uint Node1, uint Node2) const;

private:
	unsigned UnrootFromFile();
	NEWICK_TOKEN_TYPE GetTokenVerbose(TextFile &File, char szToken[],
	  unsigned uBytes) const
		{
		NEWICK_TOKEN_TYPE NTT = GetToken(File, szToken, uBytes);
		Log("GetToken %10.10s  %s\n", NTTStr(NTT), szToken);
		return NTT;
		}

	void InitCache(unsigned uCacheCount);
	void ExpandCache();
	NEWICK_TOKEN_TYPE GetToken(TextFile &File, char szToken[], unsigned uBytes) const;
	bool GetGroupFromFile(TextFile &File, unsigned uNodeIndex, double *ptrdEdgeLength);
	unsigned GetLeafCountUnrooted(unsigned uNodeIndex1, unsigned uNodeIndex2,
	  double *ptrdTotalDistance) const;
	void ToFileNodeRooted(TextFile &File, unsigned uNodeIndex) const;
	void ToFileNodeUnrooted(TextFile &File, unsigned uNodeIndex, unsigned uParent) const;
	void OrientParent(unsigned uNodeIndex, unsigned uParentNodeIndex);
	double FromClustNode(const Clust &C, unsigned uClustNodeIndex, unsigned uPhyNodeIndex);
	unsigned GetAnyNonLeafNode() const;

// Yuck. Data is made public for the convenience of Tree::Copy.
// There has to be a better way.
public:
	unsigned m_uNodeCount;
	unsigned m_uCacheCount;
	unsigned *m_uNeighbor1;
	unsigned *m_uNeighbor2;
	unsigned *m_uNeighbor3;
	double *m_dEdgeLength1;
	double *m_dEdgeLength2;
	double *m_dEdgeLength3;
	double *m_dHeight;
	bool *m_bHasEdgeLength1;
	bool *m_bHasEdgeLength2;
	bool *m_bHasEdgeLength3;
	bool *m_bHasHeight;
	unsigned *m_Ids;
	char **m_ptrName;
	bool m_bRooted;
	unsigned m_uRootNodeIndex;
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
	unsigned m_uNodeIndex1;
	unsigned m_uNodeIndex2;
	};

const unsigned NODE_CHANGED = (unsigned) (~0);

#endif // tree_h
