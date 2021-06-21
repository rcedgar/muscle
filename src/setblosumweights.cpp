/***
Code for implementing HMMer's "BLOSUM weighting" algorithm.

The algorithm was deduced by reverse-engineering the HMMer code.

The HMMer documentation refers to BLOSUM weighting as "Henikoff
simple filter weighting"

The name BLOSUM implied to me that HMMer would be using a
substitution probability matrix to compute distances, but this
turned out not to be the case.

It is notable, not to say puzzling, that the HMMer BLOSUM weighting
algorithm  is guaranteed to produce an integral NIC (number-of-indepdent-
counts, also known as effective sequence count). Presumably Eddy must
have known this, though he doesn't comment on it and he computes & stores
the value in a float.

Here's the algorithm:

Distances between two sequences are based on the average of a simple 
binary equal (one) / not equal (zero) at each position. The only thing
that has  anything to do with BLOSUM in this calculation is an obscure
(to me) constant  value of 0.62. The sequences are clustered using this
distance. If the pairwise identity (fraction of  identical positions)
is less than 0.62, they get assigned to disjoint clusters, the final
number of disjoint clusters is the NIC. This makes some intuitive sense:
I would interpret this by saying that if a set of sequences are close
enough they count as one sequence. The weight for each sequence within a
disjoint cluster is then determined to be 1 / (clustersize), from which it
follows that the sum of all weights is equal to the number of disjoint
clusters and is thus guaranteed to be an integer value. It is exactly this
sum that HMMer uses for the NIC, by default.

The individual BLOSUM sequence weights are not used for anything else in
HMMer, unless you specify that BLOSUM weighting should override the default
GSC  weighting. GSC weighting uses a different clustering algorithm to
determine  relative weights. The BLOSUM NIC is then distributed over the
GSC tree according to those relative weights.
***/

#include "muscle.h"
#include "msa.h"
#include "cluster.h"
#include "distfunc.h"

// Set weights of all sequences in the subtree under given node.
void MSA::SetBLOSUMSubtreeWeight(const ClusterNode *ptrNode, double dWeight) const
	{
	if (0 == ptrNode)
		return;

	const ClusterNode *ptrRight = ptrNode->GetRight();
	const ClusterNode *ptrLeft = ptrNode->GetLeft();

// If leaf, set weight
	if (0 == ptrRight && 0 == ptrLeft)
		{
		unsigned uIndex = ptrNode->GetIndex();
		WEIGHT w = DoubleToWeight(dWeight);
		m_Weights[uIndex] = w;
		return;
		}

// Otherwise, recursively set subtrees
	SetBLOSUMSubtreeWeight(ptrLeft, dWeight);
	SetBLOSUMSubtreeWeight(ptrRight, dWeight);
	}

// Traverse a subtree looking for clusters where all
// the leaves are sufficiently similar that they
// should be weighted as a group, i.e. given a weight
// of 1/N where N is the cluster size. The idea is
// to avoid sample bias where we have closely related
// sequences in the input alignment.
// The weight at a node is the distance between
// the two closest sequences in the left and right
// subtrees under that node. "Sufficiently similar"
// is defined as being where that minimum distance
// is less than the dMinDist threshhold. I don't know
// why the clustering is done using a minimum rather
// than a maximum or average, either of which would
// seem more natural to me.
// Return value is number of groups under this node.
// A "group" is the cluster found under a node with a
// weight less than the minimum.
unsigned MSA::SetBLOSUMNodeWeight(const ClusterNode *ptrNode, double dMinDist) const
	{
	if (0 == ptrNode)
		return 0;

	if (ptrNode->GetWeight() < dMinDist)
		{
		unsigned uClusterSize = ptrNode->GetClusterSize(); 
		assert(uClusterSize > 0);
		double dWeight = 1.0 / uClusterSize;
		SetBLOSUMSubtreeWeight(ptrNode, dWeight);
		return 1;
		}

	const ClusterNode *ptrLeft = ptrNode->GetLeft();
	const ClusterNode *ptrRight = ptrNode->GetRight();

	unsigned uLeftGroupCount = SetBLOSUMNodeWeight(ptrLeft, dMinDist);
	unsigned uRightGroupCount = SetBLOSUMNodeWeight(ptrRight, dMinDist);

	return uLeftGroupCount + uRightGroupCount;
	}

// Return value is the group count, i.e. the effective number
// of distinctly different sequences.
unsigned MSA::CalcBLOSUMWeights(ClusterTree &BlosumCluster) const
	{
// Build distance matrix
	DistFunc DF;
	unsigned uSeqCount = GetSeqCount();
	DF.SetCount(uSeqCount);
	for (unsigned i = 0; i < uSeqCount; ++i)
		for (unsigned j = i+1; j < uSeqCount; ++j)
			{
			double dDist = GetPctIdentityPair(i, j);
			assert(dDist >= 0.0 && dDist <= 1.0);
			DF.SetDist(i, j, (float) (1.0 - dDist));
			}

// Cluster based on the distance function
	BlosumCluster.Create(DF);

// Return value is HMMer's "effective sequence count".
	return SetBLOSUMNodeWeight(BlosumCluster.GetRoot(), 1.0 - BLOSUM_DIST);
	}
