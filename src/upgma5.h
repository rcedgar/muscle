#pragma once

class UPGMA5
	{
public:
	uint m_LeafCount = 0;
	uint m_TriangleSize = 0;
	uint m_InternalNodeCount = 0;
	uint m_InternalNodeIndex = 0;

// Triangular distance matrix is m_Dist, which is allocated
// as a one-dimensional vector of length m_TriangleSize.
// TriangleSubscript(i,j) maps row,column=i,j to the subscript
// into this vector.
// Row / column coordinates are a bit messy.
// Initially they are leaf indexes 0..N-1.
// But each time we create a new node (=new cluster, new subtree),
// we re-use one of the two rows that become available (the children
// of the new node). This saves memory.
// We keep track of this through the m_NodeIndex vector.
	float *m_Dist = 0;

// Distance to nearest neighbor in row i of distance matrix.
// Subscript is distance matrix row.
	float *m_MinDist = 0;

	// Nearest neighbor to row i of distance matrix.
	// Subscript is distance matrix row.
	uint *m_NearestNeighbor = 0;

// Node index of row i in distance matrix.
// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
// Subscript is distance matrix row.
	uint *m_NodeIndex = 0;

// The following vectors are defined on internal nodes,
// subscripts are internal node index 0..N-2.
// For m_Left/Right, value is the node index 0 .. 2N-2
// because a child can be internal or leaf.
	uint *m_Left = 0;
	uint *m_Right = 0;
	float *m_Height = 0;
	float *m_LeftLength = 0;
	float *m_RightLength = 0;

	vector<string> m_Labels;
	vector<vector<float> > m_DistMx;
	map<string, uint> m_LabelToIndex;

public:
	void Clear();
	void Init(const vector<string> &Labels,
	  const vector<vector<float> > &DistMx);
	void Run(const string &sLinkage, Tree &tree);
	void Run(LINKAGE Linkage, Tree &tree);
	void ReadDistMx(const string &FileName);
	void ReadDistMx2(const string &FileName);
	void ScaleDistMx(bool InputIsSimilarity = true);
	void FixEADistMx();
	void LogMe() const;
	void AddLabel(const string &Label);
	uint GetLabelIndex(const string &Label) const;

	uint TriangleSubscript(uint uIndex1, uint uIndex2) const
		{
		uint v;
		if (uIndex1 >= uIndex2)
			v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
		else
			v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
		assert(v < (m_LeafCount*(m_LeafCount - 1))/2);
		return v;
		}
	};
