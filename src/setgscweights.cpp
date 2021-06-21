/***
Gerstein/Sonnhammer/Chothia ad hoc sequence weighting.
The algorithm was deduced by reverse-engineering the
HMMer code.

I used an alternative representation that I prefer over
HMMer's. The HMMer code is full of tree manipulations
that do something to the left child and then the equivalent
thing to the right child. It was clear that there must be
a re-formulation that does everything once for each node,
which would reduce the number of operations expressed
in the code by a factor of two. This gives a more elegant
and less error-prone way to code it.

These notes explain the correspondence between my design
and Eddy's.

HMMer stores a data structure phylo_s for each non-leaf
node in the cluster tree. This structure contains the
following fields:

	diff		Weight of the node
	lblen		Left branch length
	rblen		Right branch length

The lblen and rblen branch lengths are calculated as:

	this.lblen = this.diff - left.diff
	this.rblen = this.diff - right.diff

My code stores one ClusterNode data structure per node
in the cluster tree, including leaves. I store only the
weight. I can recover the HMMer branch length fields
in a trivial O(1) calculation as follows:

	lblen = Node.GetWeight() - Node.GetLeft()->GetWeight()
	rblen = Node.GetWeight() - Node.GetRight()->GetWeight()

For the GSC weights calculation, HMMer constructs the
following vectors, which have entries for all nodes,
including leaves:

	lwt		Left weight
	rwt		Right weight

The "left weight" is calculated as the sum of the weights in
all the nodes reachable through the left branch, including
the node itself. (This is not immediately obvious from the
code, which does the calculation using branch lengths rather
than weights, but this is an equivalent, and to my mind clearer,
statement of what they are). Similarly, the "right weight" is
the sum of all weights reachable via the right branch. I define
the "cluster weight" to be the summed weight of all nodes in the
subtree under the node, including the node itself. I provide
a function Node.GetClusterWeight() which calculates the cluster
weight using a O(ln N) recursion through the tree. The lwt and
rwt values can be recovered as follows:

	lwt		= Node.GetLeft()->GetClusterWeight()
			+ Node.GetWeight()

	lwt		= Node.GetLeft()->GetClusterWeight()
			+ Node.GetWeight()

HMMer calculates a further vector fwt as follows.

	this.fwt = parent.fwt * parent.lwt / (parent.lwt + parent.rwt)

This applies to nodes reached via a left branch, for nodes reached
via a right branch:

	this.fwt = parent.fwt * parent.rwt / (parent.lwt + parent.rwt)

The values of fwt at the leaf nodes are the final GSC weights.
We derive the various terms using our equivalents.

	parent.lwt	= Parent.GetLeft()->GetClusterWeight()
				+ Parent.GetWeight()

	parent.rwt	= Parent.GetRight()->GetClusterWeight()
				+ Parent.GetWeight()

	parent.lwt + parent.rwt =
				{ Parent.GetLeft()->GetClusterWeight()
				+ Parent.GetRight()->GetClusterWeight()
				+ Parent.GetWeight() }
				+ Parent.GetWeight()

We recognize the term {...} as the cluster weight of the
parent, so

	parent.lwt + parent.rwt
				= Parent.GetClusterWeight()
				+ Parent.GetWeight()

As you would expect, repeating this exercise for parent.rwt gives
exactly the same expression.

The GSC weights (fwt) are stored in the Weight2 field of the cluster
tree, the Weight field stores the original (BLOSUM) weights used
as input to this algorithm.
***/

#include "muscle.h"
#include "msa.h"
#include "cluster.h"
#include "distfunc.h"

// Set weights of all sequences in the subtree under given node.
void MSA::SetSubtreeWeight2(const ClusterNode *ptrNode) const
	{
	if (0 == ptrNode)
		return;

	const ClusterNode *ptrRight = ptrNode->GetRight();
	const ClusterNode *ptrLeft = ptrNode->GetLeft();

// If leaf, set weight
	if (0 == ptrRight && 0 == ptrLeft)
		{
		unsigned uIndex = ptrNode->GetIndex();
		double dWeight = ptrNode->GetWeight2();
		WEIGHT w = DoubleToWeight(dWeight);
		m_Weights[uIndex] = w;
		return;
		}

// Otherwise, recursively set subtrees
	SetSubtreeWeight2(ptrLeft);
	SetSubtreeWeight2(ptrRight);
	}

void MSA::SetSubtreeGSCWeight(ClusterNode *ptrNode) const
	{
	if (0 == ptrNode)
		return;

	ClusterNode *ptrParent = ptrNode->GetParent();
	double dParentWeight2 = ptrParent->GetWeight2();
	double dParentClusterWeight = ptrParent->GetClusterWeight();
	if (0.0 == dParentClusterWeight)
		{
		double dThisClusterSize = ptrNode->GetClusterSize();
		double dParentClusterSize = ptrParent->GetClusterSize();
		double dWeight2 =
		  dParentWeight2*dThisClusterSize/dParentClusterSize;
		ptrNode->SetWeight2(dWeight2);
		}
	else
		{
	// Could cache cluster weights for better performance.
	// We calculate cluster weight of each node twice, so this
	// would give x2 improvement.
	// As weighting is not very expensive, we don't care.
		double dThisClusterWeight = ptrNode->GetClusterWeight();
		double dParentWeight = ptrParent->GetWeight();

		double dNum = dThisClusterWeight + dParentWeight;
		double dDenom = dParentClusterWeight + dParentWeight;
		double dWeight2 = dParentWeight2*(dNum/dDenom);

		ptrNode->SetWeight2(dWeight2);
		}

	SetSubtreeGSCWeight(ptrNode->GetLeft());
	SetSubtreeGSCWeight(ptrNode->GetRight());
	}

void MSA::SetGSCWeights() const
	{
	ClusterTree CT;
	CalcBLOSUMWeights(CT);

// Calculate weights and store in tree.
	ClusterNode *ptrRoot = CT.GetRoot();
	ptrRoot->SetWeight2(1.0);
	SetSubtreeGSCWeight(ptrRoot->GetLeft());
	SetSubtreeGSCWeight(ptrRoot->GetRight());

// Copy weights from tree to MSA.
	SetSubtreeWeight2(ptrRoot);
	}
 
void MSA::ListWeights() const
	{
	const unsigned uSeqCount = GetSeqCount();
	Log("Weights:\n");
	WEIGHT wTotal = 0;
	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		wTotal += GetSeqWeight(n);
		Log("%6.3f %s\n", GetSeqWeight(n), GetSeqName(n));
		}
	Log("Total weights = %6.3f, should be 1.0\n", wTotal);
	}
