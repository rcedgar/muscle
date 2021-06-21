#ifndef ClustSetDF_h
#define ClustSetDF_h

class MSA;
class Clust;

#include "clustset.h"
#include "distfunc.h"
#include "msa.h"

class ClustSetDF : public ClustSet
	{
public:
	ClustSetDF(const DistFunc &DF) :
		m_ptrDF(&DF)
		{
		}

public:
	virtual unsigned GetLeafCount()
		{
		return m_ptrDF->GetCount();
		}
	virtual const char *GetLeafName(unsigned uNodeIndex)
		{
		return m_ptrDF->GetName(uNodeIndex);
		}
	virtual unsigned GetLeafId(unsigned uNodeIndex)
		{
		return m_ptrDF->GetId(uNodeIndex);
		}
	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
	  double *ptrdLeftLength, double *ptrdRightLength)
		{
		Quit("ClustSetDF::JoinNodes, should never be called");
		}
	virtual double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2)
		{
		return m_ptrDF->GetDist(uNodeIndex1, uNodeIndex2);
		}

private:
	const DistFunc *m_ptrDF;
	};

#endif	// ClustSetDF_h
