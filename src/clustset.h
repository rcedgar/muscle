#ifndef ClustSet_h
#define ClustSet_h

enum JOIN;
enum LINKAGE;
class Clust;

class ClustSet
	{
public:
	virtual unsigned GetLeafCount() = 0;
	virtual double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2) = 0;
	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
	  double *ptrdLeftLength, double *ptrdRightLength) = 0;
	virtual const char *GetLeafName(unsigned uNodeIndex) = 0;
	virtual unsigned GetLeafId(unsigned uNodeIndex) = 0;
	};

#endif	// ClustSet_h
