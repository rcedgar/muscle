#ifndef EdgeList_h
#define EdgeList_h

class EdgeList
	{
public:
	EdgeList();
	virtual ~EdgeList();

public:
	void Clear();
	void Add(unsigned uNode1, unsigned uNode2);
	unsigned GetCount() const;
	void GetEdge(unsigned uIndex, unsigned *ptruNode1, unsigned *ptruNode2) const;
	void Copy(const EdgeList &rhs);
	void LogMe() const;

private:
	void Expand();

private:
	unsigned m_uCount;
	unsigned m_uCacheSize;
	unsigned *m_uNode1;
	unsigned *m_uNode2;
	};

#endif	// EdgeList_h
