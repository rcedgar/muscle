#include "muscle.h"
#include "tree.h"
#include "seqvect.h"
#include "profile.h"
#include "msa.h"
#include "pwpath.h"
#include "estring.h"

#define TRACE		0
#define VALIDATE	0

static void PathSeq(const Seq &s, const PWPath &Path, bool bRight, Seq &sOut)
	{
	int *esA;
	int *esB;
	PathToEstrings(Path, &esA, &esB);

	const unsigned uSeqLength = s.Length();
	const unsigned uEdgeCount = Path.GetEdgeCount();

	sOut.Clear();
	sOut.SetName(s.GetName());
	unsigned uPos = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		char cType = Edge.cType;
		if (bRight)
			{
			if (cType == 'I')
				cType = 'D';
			else if (cType == 'D')
				cType = 'I';
			}
		switch (cType)
			{
		case 'M':
			sOut.AppendChar(s[uPos++]);
			break;
		case 'D':
			sOut.AppendChar('-');
			break;
		case 'I':
			sOut.AppendChar(s[uPos++]);
			break;
		default:
			Quit("PathSeq, invalid edge type %c", cType);
			}
		}
	}

#if	VALIDATE

static void MakeRootSeq(const Seq &s, const Tree &GuideTree, unsigned uLeafNodeIndex,
  const ProgNode Nodes[], Seq &sRoot)
	{
	sRoot.Copy(s);
	unsigned uNodeIndex = uLeafNodeIndex;
	for (;;)
		{
	  	unsigned uParent = GuideTree.GetParent(uNodeIndex);
		if (NULL_NEIGHBOR == uParent)
			break;
		bool bRight = (GuideTree.GetLeft(uParent) == uNodeIndex);
		uNodeIndex = uParent;
		const PWPath &Path = Nodes[uNodeIndex].m_Path;
		Seq sTmp;
		PathSeq(sRoot, Path, bRight, sTmp);
		sTmp.SetId(0);
		sRoot.Copy(sTmp);
		}
	}

#endif	// VALIDATE

static int *MakeRootSeqE(const Seq &s, const Tree &GuideTree, unsigned uLeafNodeIndex,
  const ProgNode Nodes[], Seq &sRoot, int *Estring1, int *Estring2)
	{
	int *EstringCurr = Estring1;
	int *EstringNext = Estring2;

	const unsigned uSeqLength = s.Length();
	EstringCurr[0] = uSeqLength;
	EstringCurr[1] = 0;

	unsigned uNodeIndex = uLeafNodeIndex;
	for (;;)
		{
	  	unsigned uParent = GuideTree.GetParent(uNodeIndex);
		if (NULL_NEIGHBOR == uParent)
			break;
		bool bRight = (GuideTree.GetLeft(uParent) == uNodeIndex);
		uNodeIndex = uParent;
		const PWPath &Path = Nodes[uNodeIndex].m_Path;
		const int *EstringNode = bRight ?
		  Nodes[uNodeIndex].m_EstringL : Nodes[uNodeIndex].m_EstringR;

		MulEstrings(EstringCurr, EstringNode, EstringNext);
#if	TRACE
		Log("\n");
		Log("Curr=");
		LogEstring(EstringCurr);
		Log("\n");
		Log("Node=");
		LogEstring(EstringNode);
		Log("\n");
		Log("Prod=");
		LogEstring(EstringNext);
		Log("\n");
#endif
		int *EstringTmp = EstringNext;
		EstringNext = EstringCurr;
		EstringCurr = EstringTmp;
		}
	EstringOp(EstringCurr, s, sRoot);

#if	TRACE
	Log("Root estring=");
	LogEstring(EstringCurr);
	Log("\n");
	Log("Root seq=");
	sRoot.LogMe();
#endif
	return EstringCurr;
	}

static unsigned GetFirstNodeIndex(const Tree &tree)
	{
	if (g_bStable)
		return 0;
	return tree.FirstDepthFirstNode();
	}

static unsigned GetNextNodeIndex(const Tree &tree, unsigned uPrevNodeIndex)
	{
	if (g_bStable)
		{
		const unsigned uNodeCount = tree.GetNodeCount();
		unsigned uNodeIndex = uPrevNodeIndex;
		for (;;)
			{
			++uNodeIndex;
			if (uNodeIndex >= uNodeCount)
				return NULL_NEIGHBOR;
			if (tree.IsLeaf(uNodeIndex))
				return uNodeIndex;
			}
		}
	unsigned uNodeIndex = uPrevNodeIndex;
	for (;;)
		{
		uNodeIndex = tree.NextDepthFirstNode(uNodeIndex);
		if (NULL_NEIGHBOR == uNodeIndex || tree.IsLeaf(uNodeIndex))
			return uNodeIndex;
		}
	}

void MakeRootMSA(const SeqVect &v, const Tree &GuideTree, ProgNode Nodes[],
  MSA &a)
	{
#if	TRACE
	Log("MakeRootMSA Tree=");
	GuideTree.LogMe();
#endif
	const unsigned uSeqCount = v.GetSeqCount();
	unsigned uColCount = uInsane;
	unsigned uSeqIndex = 0;
	const unsigned uTreeNodeCount = GuideTree.GetNodeCount();
	const unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	const PWPath &RootPath = Nodes[uRootNodeIndex].m_Path;
	const unsigned uRootColCount = RootPath.GetEdgeCount();
	const unsigned uEstringSize = uRootColCount + 1;
	int *Estring1 = new int[uEstringSize];
	int *Estring2 = new int[uEstringSize];
	SetProgressDesc("Root alignment");

	unsigned uTreeNodeIndex = GetFirstNodeIndex(GuideTree);
	do
		{
		Progress(uSeqIndex, uSeqCount);

		unsigned uId = GuideTree.GetLeafId(uTreeNodeIndex);
		const Seq &s = *(v[uId]);

		Seq sRootE;
		int *es = MakeRootSeqE(s, GuideTree, uTreeNodeIndex, Nodes, sRootE,
		  Estring1, Estring2);
		Nodes[uTreeNodeIndex].m_EstringL = EstringNewCopy(es);

#if	VALIDATE
		Seq sRoot;
		MakeRootSeq(s, GuideTree, uTreeNodeIndex, Nodes, sRoot);
		if (!sRoot.Eq(sRootE))
			{
			Log("sRoot=");
			sRoot.LogMe();
			Log("sRootE=");
			sRootE.LogMe();
			Quit("Root seqs differ");
			}
#if	TRACE
		Log("MakeRootSeq=\n");
		sRoot.LogMe();
#endif
#endif

		if (uInsane == uColCount)
			{
			uColCount = sRootE.Length();
			a.SetSize(uSeqCount, uColCount);
			}
		else
			{
			assert(uColCount == sRootE.Length());
			}
		a.SetSeqName(uSeqIndex, s.GetName());
		a.SetSeqId(uSeqIndex, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			a.SetChar(uSeqIndex, uColIndex, sRootE[uColIndex]);
		++uSeqIndex;

		uTreeNodeIndex = GetNextNodeIndex(GuideTree, uTreeNodeIndex);
		}
	while (NULL_NEIGHBOR != uTreeNodeIndex);

	delete[] Estring1;
	delete[] Estring2;

	ProgressStepsDone();
	assert(uSeqIndex == uSeqCount);
	}
