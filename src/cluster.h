class DistFunc;

class ClusterNode
	{
	friend class ClusterTree;
public:
	ClusterNode()
		{
		m_dWeight = 0.0;
		m_dWeight2 = 0.0;
		m_ptrLeft = 0;
		m_ptrRight = 0;
		m_ptrParent = 0;
		m_uIndex = 0;
		m_ptrPrevDisjoint = 0;
		m_ptrNextDisjoint = 0;
		}
	~ClusterNode() {}

public:
	unsigned GetIndex() const { return m_uIndex; }
	ClusterNode *GetLeft() const { return m_ptrLeft; }
	ClusterNode *GetRight() const { return m_ptrRight; }
	ClusterNode *GetParent() const { return m_ptrParent; }
	double GetWeight() const { return m_dWeight; }

	const ClusterNode *GetClusterLeaf(unsigned uLeafIndex) const;
	unsigned GetClusterSize() const;
	double GetClusterWeight() const;
	double GetLeftBranchWeight() const;
	double GetRightBranchWeight() const;
	double GetLeftWeight() const;
	double GetRightWeight() const;

	void LogMe() const;

	double GetWeight2() const { return m_dWeight2; }
	void SetWeight2(double dWeight2) { m_dWeight2 = dWeight2; }

protected:
	void SetIndex(unsigned uIndex) { m_uIndex = uIndex; }
	void SetWeight(double dWeight) { m_dWeight = dWeight; }
	void SetLeft(ClusterNode *ptrLeft) { m_ptrLeft = ptrLeft; }
	void SetRight(ClusterNode *ptrRight) { m_ptrRight = ptrRight; }
	void SetParent(ClusterNode *ptrParent) { m_ptrParent = ptrParent; }
	void SetNextDisjoint(ClusterNode *ptrNode) { m_ptrNextDisjoint = ptrNode; }
	void SetPrevDisjoint(ClusterNode *ptrNode) { m_ptrPrevDisjoint = ptrNode; }

	ClusterNode *GetNextDisjoint() { return m_ptrNextDisjoint; }
	ClusterNode *GetPrevDisjoint() { return m_ptrPrevDisjoint; }

private:
	double m_dWeight;
	double m_dWeight2;
	unsigned m_uIndex;
	ClusterNode *m_ptrLeft;
	ClusterNode *m_ptrRight;
	ClusterNode *m_ptrParent;
	ClusterNode *m_ptrNextDisjoint;
	ClusterNode *m_ptrPrevDisjoint;
	};

class ClusterTree
	{
public:
	ClusterTree();
	virtual ~ClusterTree();

	void Create(const DistFunc &DF);

	ClusterNode *GetRoot() const;
	void LogMe() const;

protected:
	void Join(ClusterNode *ptrNode1, ClusterNode *ptrNode2,
	  ClusterNode *ptrJoin);
	void AddToDisjoints(ClusterNode *ptrNode);
	void DeleteFromDisjoints(ClusterNode *ptrNode);
	void Validate(unsigned uNodeCount);

private:
	ClusterNode *m_ptrDisjoints;
	ClusterNode *m_Nodes;
	unsigned m_uNodeCount;
	unsigned m_uLeafCount;
	};
