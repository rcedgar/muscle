#ifndef Clust_h
#define Clust_h

class Clust;
class ClustNode;
class ClustSet;
class Phylip;
class SortedNode;

const unsigned RB_NIL = ((unsigned) 0xfff0);

class ClustNode
	{
public:
	ClustNode()
		{
		m_uIndex = uInsane;
		m_uSize = uInsane;
		m_dLength = (float) dInsane;
		m_ptrLeft = 0;
		m_ptrRight = 0;
		m_ptrParent = 0;
		m_ptrNextCluster = 0;
		m_ptrPrevCluster = 0;
		m_uLeafIndexes = 0;
		}
	~ClustNode()
		{
		delete[] m_uLeafIndexes;
		}
	unsigned m_uIndex;
	unsigned m_uSize;
	float m_dLength;
	ClustNode *m_ptrLeft;
	ClustNode *m_ptrRight;
	ClustNode *m_ptrParent;
	ClustNode *m_ptrNextCluster;
	ClustNode *m_ptrPrevCluster;
	unsigned *m_uLeafIndexes;
	};

class Clust
	{
public:
	Clust();
	virtual ~Clust();

	void Create(ClustSet &Set, CLUSTER Method);

	unsigned GetLeafCount() const;

	unsigned GetClusterCount() const;
	unsigned GetClusterSize(unsigned uNodeIndex) const;
	unsigned GetLeaf(unsigned uClusterIndex, unsigned uLeafIndex) const;

	unsigned GetNodeCount() const { return 2*m_uLeafCount - 1; }
	const ClustNode &GetRoot() const { return m_Nodes[GetRootNodeIndex()]; }
	unsigned GetRootNodeIndex() const { return m_uNodeCount - 1; }

	const ClustNode &GetNode(unsigned uNodeIndex) const;
	bool IsLeaf(unsigned uNodeIndex) const;
	unsigned GetLeftIndex(unsigned uNodeIndex) const;
	unsigned GetRightIndex(unsigned uNodeIndex) const;
	float GetLength(unsigned uNodeIndex) const;
	float GetHeight(unsigned uNodeIndex) const;
	const char *GetNodeName(unsigned uNodeIndex) const;
	unsigned GetNodeId(unsigned uNodeIndex) const;

	JOIN GetJoinStyle() const { return m_JoinStyle; }
	LINKAGE GetCentroidStyle() const { return m_CentroidStyle; }

	void SetDist(unsigned uIndex1, unsigned uIndex2, float dDist);
	float GetDist(unsigned uIndex1, unsigned uIndex2) const;

	void ToPhylip(Phylip &tree);

	void LogMe() const;

//private:
	void SetLeafCount(unsigned uLeafCount);

	void CreateCluster();
	void JoinNodes(unsigned uLeftNodeIndex, unsigned uRightNodeIndex, 
	  float dLeftLength, float dRightLength, unsigned uNewNodeIndex);

	void ChooseJoin(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  float *ptrdLeftLength, float *ptrdRightLength);
	void ChooseJoinNeighborJoining(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  float *ptrdLeftLength, float *ptrdRightLength);
	void ChooseJoinNearestNeighbor(unsigned *ptruLeftIndex, unsigned *ptruRightIndex,
	  float *ptrdLeftLength, float *ptrdRightLength);

	float ComputeDist(unsigned uNewNodeIndex, unsigned uNodeIndex);
	float ComputeDistAverageLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	float ComputeDistMinLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	float ComputeDistMaxLinkage(unsigned uNewNodeIndex, unsigned uNodeIndex);
	float ComputeDistNeighborJoining(unsigned uNewNewIndex, unsigned uNodeIndex);
	float ComputeDistMAFFT(unsigned uNewNewIndex, unsigned uNodeIndex);

	float Calc_r(unsigned uNodeIndex) const;

	unsigned VectorIndex(unsigned uIndex1, unsigned uIndex2) const;

	unsigned GetFirstCluster() const;
	unsigned GetNextCluster(unsigned uNodeIndex) const;

	float ComputeMetric(unsigned uIndex1, unsigned uIndex2) const;
	float ComputeMetricNearestNeighbor(unsigned i, unsigned j) const;
	float ComputeMetricNeighborJoining(unsigned i, unsigned j) const;

	void InitMetric(unsigned uMaxNodeIndex);
	void InsertMetric(unsigned uIndex1, unsigned uIndex2, float dMetric);
	float GetMinMetric(unsigned *ptruIndex1, unsigned *ptruIndex2) const;
	float GetMinMetricBruteForce(unsigned *ptruIndex1, unsigned *ptruIndex2) const;
	void DeleteMetric(unsigned uIndex);
	void DeleteMetric(unsigned uIndex1, unsigned uIndex2);
	void ListMetric() const;

	void DeleteFromClusterList(unsigned uNodeIndex);
	void AddToClusterList(unsigned uNodeIndex);

	void RBDelete(unsigned RBNode);
	unsigned RBInsert(unsigned i, unsigned j, float fMetric);

	unsigned RBNext(unsigned RBNode) const;
	unsigned RBPrev(unsigned RBNode) const;
	unsigned RBMin(unsigned RBNode) const;
	unsigned RBMax(unsigned RBNode) const;

	void ValidateRB(const char szMsg[] = 0) const;
	void ValidateRBNode(unsigned Node, const char szMsg[]) const;

//private:
	JOIN m_JoinStyle;
	LINKAGE m_CentroidStyle;
	ClustNode *m_Nodes;
	unsigned *m_ClusterIndexToNodeIndex;
	unsigned *m_NodeIndexToClusterIndex;
	unsigned m_uLeafCount;
	unsigned m_uNodeCount;
	unsigned m_uClusterCount;
	unsigned m_uTriangularMatrixSize;
	float *m_dDist;
	ClustSet *m_ptrSet;
	ClustNode *m_ptrClusterList;
	};

#endif // Clust_h
