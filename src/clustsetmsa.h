#ifndef ClustSetMSA_h
#define ClustSetMSA_h

class MSA;
class Clust;

#include "clustset.h"
#include "msadist.h"

// Distance matrix based set.
// Computes distances between leaves, never between
// joined clusters (leaves this to distance matrix method).
class ClustSetMSA : public ClustSet
	{
public:
	ClustSetMSA(const MSA &msa, MSADist &MD) :
		m_ptrMSA(&msa),
		m_ptrMSADist(&MD)
		{
		}

public:
	virtual unsigned GetLeafCount()
		{
		return m_ptrMSA->GetSeqCount();
		}
	virtual const char *GetLeafName(unsigned uNodeIndex)
		{
		return m_ptrMSA->GetSeqName(uNodeIndex);
		}
	virtual unsigned GetLeafId(unsigned uNodeIndex)
		{
		return m_ptrMSA->GetSeqId(uNodeIndex);
		}
	virtual void JoinNodes(const Clust &C, unsigned uLeftNodeIndex,
	  unsigned uRightNodeIndex, unsigned uJoinedNodeIndex,
	  double *ptrdLeftLength, double *ptrdRightLength)
		{
		Quit("ClustSetMSA::JoinNodes, should never be called");
		}
	virtual double ComputeDist(const Clust &C, unsigned uNodeIndex1,
	  unsigned uNodeIndex2)
		{
		return m_ptrMSADist->ComputeDist(*m_ptrMSA, uNodeIndex1, uNodeIndex2);
		}

public:
	const MSA &GetMSA();

private:
	const MSA *m_ptrMSA;
	MSADist *m_ptrMSADist;
	};

#endif	// ClustSetMSA_h
