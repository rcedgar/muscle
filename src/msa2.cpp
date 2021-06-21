#include "muscle.h"
#include "msa.h"
#include "seqvect.h"
#include "profile.h"
#include "tree.h"

// These global variables are a hack to allow the tree
// dependent iteration code to communicate the edge
// used to divide the tree. The three-way weighting
// scheme needs to know this edge in order to compute
// sequence weights.
static const Tree *g_ptrMuscleTree = 0;
unsigned g_uTreeSplitNode1 = NULL_NEIGHBOR;
unsigned g_uTreeSplitNode2 = NULL_NEIGHBOR;

void MSA::GetFractionalWeightedCounts(unsigned uColIndex, bool bNormalize,
  FCOUNT fcCounts[], FCOUNT *ptrfcGapStart, FCOUNT *ptrfcGapEnd,
  FCOUNT *ptrfcGapExtend, FCOUNT *ptrfOcc,
  FCOUNT *ptrfcLL, FCOUNT *ptrfcLG, FCOUNT *ptrfcGL, FCOUNT *ptrfcGG) const
	{
	const unsigned uSeqCount = GetSeqCount();
	const unsigned uColCount = GetColCount();

	memset(fcCounts, 0, g_AlphaSize*sizeof(FCOUNT));
	WEIGHT wTotal = 0;
	FCOUNT fGap = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const WEIGHT w = GetSeqWeight(uSeqIndex);
		if (IsGap(uSeqIndex, uColIndex))
			{
			fGap += w;
			continue;
			}
		else if (IsWildcard(uSeqIndex, uColIndex))
			{
			const unsigned uLetter = GetLetterEx(uSeqIndex, uColIndex);
			switch (g_Alpha)
				{
			case ALPHA_Amino:
				switch (uLetter)
					{
				case AX_B:		// D or N
					fcCounts[AX_D] += w/2;
					fcCounts[AX_N] += w/2;
					break;
				case AX_Z:		// E or Q
					fcCounts[AX_E] += w/2;
					fcCounts[AX_Q] += w/2;
					break;
				default:		// any
					{
					const FCOUNT f = w/20;
					for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
						fcCounts[uLetter] += f;
					break;
					}
					}
				break;

			case ALPHA_DNA:
			case ALPHA_RNA:
				switch (uLetter)
					{
				case AX_R:	// G or A
					fcCounts[NX_G] += w/2;
					fcCounts[NX_A] += w/2;
					break;
				case AX_Y:	// C or T/U
					fcCounts[NX_C] += w/2;
					fcCounts[NX_T] += w/2;
					break;
				default:	// any
					const FCOUNT f = w/20;
					for (unsigned uLetter = 0; uLetter < 4; ++uLetter)
						fcCounts[uLetter] += f;
					break;
					}
				break;

			default:
				Quit("Alphabet %d not supported", g_Alpha);
				}
			continue;
			}
		unsigned uLetter = GetLetter(uSeqIndex, uColIndex);
		fcCounts[uLetter] += w;
		wTotal += w;
		}
	*ptrfOcc = (float) (1.0 - fGap);

	if (bNormalize && wTotal > 0)
		{
		if (wTotal > 1.001)
			Quit("wTotal=%g\n", wTotal);
		for (unsigned uLetter = 0; uLetter < g_AlphaSize; ++uLetter)
			fcCounts[uLetter] /= wTotal;
//		AssertNormalized(fcCounts);
		}

	FCOUNT fcStartCount = 0;
	if (uColIndex == 0)
		{
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			if (IsGap(uSeqIndex, uColIndex))
				fcStartCount += GetSeqWeight(uSeqIndex);
		}
	else
		{
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			if (IsGap(uSeqIndex, uColIndex) && !IsGap(uSeqIndex, uColIndex - 1))
				fcStartCount += GetSeqWeight(uSeqIndex);
		}

	FCOUNT fcEndCount = 0;
	if (uColCount - 1 == uColIndex)
		{
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			if (IsGap(uSeqIndex, uColIndex))
				fcEndCount += GetSeqWeight(uSeqIndex);
		}
	else
		{
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			if (IsGap(uSeqIndex, uColIndex) && !IsGap(uSeqIndex, uColIndex + 1))
				fcEndCount += GetSeqWeight(uSeqIndex);
		}

	FCOUNT LL = 0;
	FCOUNT LG = 0;
	FCOUNT GL = 0;
	FCOUNT GG = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		WEIGHT w = GetSeqWeight(uSeqIndex);
		bool bLetterHere = !IsGap(uSeqIndex, uColIndex);
		bool bLetterPrev = (uColIndex == 0 || !IsGap(uSeqIndex, uColIndex - 1));
		if (bLetterHere)
			{
			if (bLetterPrev)
				LL += w;
			else
				GL += w;
			}
		else
			{
			if (bLetterPrev)
				LG += w;
			else
				GG += w;
			}
		}

	FCOUNT fcExtendCount = 0;
	if (uColIndex > 0 && uColIndex < GetColCount() - 1)
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			if (IsGap(uSeqIndex, uColIndex) && IsGap(uSeqIndex, uColIndex - 1) &&
			  IsGap(uSeqIndex, uColIndex + 1))
				fcExtendCount += GetSeqWeight(uSeqIndex);

	*ptrfcLL = LL;
	*ptrfcLG = LG;
	*ptrfcGL = GL;
	*ptrfcGG = GG;
	*ptrfcGapStart = fcStartCount;
	*ptrfcGapEnd = fcEndCount;
	*ptrfcGapExtend = fcExtendCount;
	}

// Return true if the given column has no gaps and all
// its residues are in the same biochemical group.
bool MSAColIsConservative(const MSA &msa, unsigned uColIndex)
	{
	extern unsigned ResidueGroup[];

	const unsigned uSeqCount = msa.GetColCount();
	if (0 == uSeqCount)
		Quit("MSAColIsConservative: empty alignment");

	if (msa.IsGap(0, uColIndex))
		return false;

	unsigned uLetter = msa.GetLetterEx(0, uColIndex);
	const unsigned uGroup = ResidueGroup[uLetter];

	for (unsigned uSeqIndex = 1; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		if (msa.IsGap(uSeqIndex, uColIndex))
			return false;
		uLetter = msa.GetLetter(uSeqIndex, uColIndex);
		if (ResidueGroup[uLetter] != uGroup)
			return false;
		}
	return true;
	}

void MSAFromSeqRange(const MSA &msaIn, unsigned uFromSeqIndex, unsigned uSeqCount,
  MSA &msaOut)
	{
	const unsigned uColCount = msaIn.GetColCount();
	msaOut.SetSize(uSeqCount, uColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const char *ptrName = msaIn.GetSeqName(uFromSeqIndex + uSeqIndex);
		msaOut.SetSeqName(uSeqIndex, ptrName);

		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uFromSeqIndex + uSeqIndex, uColIndex);
			msaOut.SetChar(uSeqIndex, uColIndex, c);
			}
		}
	}

void MSAFromColRange(const MSA &msaIn, unsigned uFromColIndex, unsigned uColCount,
  MSA &msaOut)
	{
	const unsigned uSeqCount = msaIn.GetSeqCount();
	const unsigned uInColCount = msaIn.GetColCount();

	if (uFromColIndex + uColCount - 1 > uInColCount)
		Quit("MSAFromColRange, out of bounds");

	msaOut.SetSize(uSeqCount, uColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const char *ptrName = msaIn.GetSeqName(uSeqIndex);
		unsigned uId = msaIn.GetSeqId(uSeqIndex);
		msaOut.SetSeqName(uSeqIndex, ptrName);
		msaOut.SetSeqId(uSeqIndex, uId);

		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uSeqIndex, uFromColIndex + uColIndex);
			msaOut.SetChar(uSeqIndex, uColIndex, c);
			}
		}
	}

void SeqVectFromMSA(const MSA &msa, SeqVect &v)
	{
	v.Clear();
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq s;
		msa.GetSeq(uSeqIndex, s);

		s.StripGaps();
		//if (0 == s.Length())
		//	continue;

		const char *ptrName = msa.GetSeqName(uSeqIndex);
		s.SetName(ptrName);

		v.AppendSeq(s);
		}
	}

void DeleteGappedCols(MSA &msa)
	{
	unsigned uColIndex = 0;
	for (;;)
		{
		if (uColIndex >= msa.GetColCount())
			break;
		if (msa.IsGapColumn(uColIndex))
			msa.DeleteCol(uColIndex);
		else
			++uColIndex;
		}
	}

void MSAFromSeqSubset(const MSA &msaIn, const unsigned uSeqIndexes[], unsigned uSeqCount,
  MSA &msaOut)
	{
	const unsigned uColCount = msaIn.GetColCount();
	msaOut.SetSize(uSeqCount, uColCount);
	for (unsigned uSeqIndexOut = 0; uSeqIndexOut < uSeqCount; ++uSeqIndexOut)
		{
		unsigned uSeqIndexIn = uSeqIndexes[uSeqIndexOut];
		const char *ptrName = msaIn.GetSeqName(uSeqIndexIn);
		unsigned uId = msaIn.GetSeqId(uSeqIndexIn);
		msaOut.SetSeqName(uSeqIndexOut, ptrName);
		msaOut.SetSeqId(uSeqIndexOut, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uSeqIndexIn, uColIndex);
			msaOut.SetChar(uSeqIndexOut, uColIndex, c);
			}
		}
	}

void AssertMSAEqIgnoreCaseAndGaps(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	if (uSeqCount1 != uSeqCount2)
		Quit("Seq count differs");

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount1; ++uSeqIndex)
		{
		Seq seq1;
		msa1.GetSeq(uSeqIndex, seq1);

		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);

		Seq seq2;
		msa2.GetSeq(uSeqIndex2, seq2);

		if (!seq1.EqIgnoreCaseAndGaps(seq2))
			{
			Log("Input:\n");
			seq1.LogMe();
			Log("Output:\n");
			seq2.LogMe();
			Quit("Seq %s differ ", msa1.GetSeqName(uSeqIndex));
			}
		}
	}

void AssertMSAEq(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	if (uSeqCount1 != uSeqCount2)
		Quit("Seq count differs");

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount1; ++uSeqIndex)
		{
		Seq seq1;
		msa1.GetSeq(uSeqIndex, seq1);

		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);

		Seq seq2;
		msa2.GetSeq(uSeqIndex2, seq2);

		if (!seq1.Eq(seq2))
			{
			Log("Input:\n");
			seq1.LogMe();
			Log("Output:\n");
			seq2.LogMe();
			Quit("Seq %s differ ", msa1.GetSeqName(uSeqIndex));
			}
		}
	}

void SetMSAWeightsMuscle(MSA &msa)
	{
	SEQWEIGHT Method = GetSeqWeightMethod();
	switch (Method)
		{
	case SEQWEIGHT_None:
		msa.SetUniformWeights();
		return;

	case SEQWEIGHT_Henikoff:
		msa.SetHenikoffWeights();
		return;

	case SEQWEIGHT_HenikoffPB:
		msa.SetHenikoffWeightsPB();
		return;

	case SEQWEIGHT_GSC:
		msa.SetGSCWeights();
		return;

	case SEQWEIGHT_ClustalW:
		SetClustalWWeightsMuscle(msa);
		return;
	
	case SEQWEIGHT_ThreeWay:
		SetThreeWayWeightsMuscle(msa);
		return;
		}
	Quit("SetMSAWeightsMuscle, Invalid method=%d", Method);
	}

static WEIGHT *g_MuscleWeights;
static unsigned g_uMuscleIdCount;

WEIGHT GetMuscleSeqWeightById(unsigned uId)
	{
	if (0 == g_MuscleWeights)
		Quit("g_MuscleWeights = 0");
	if (uId >= g_uMuscleIdCount)
		Quit("GetMuscleSeqWeightById(%u): count=%u",
		  uId, g_uMuscleIdCount);

	return g_MuscleWeights[uId];
	}

void SetMuscleTree(const Tree &tree)
	{
	g_ptrMuscleTree = &tree;

	if (SEQWEIGHT_ClustalW != GetSeqWeightMethod())
		return;

	delete[] g_MuscleWeights;

	const unsigned uLeafCount = tree.GetLeafCount();
	g_uMuscleIdCount = uLeafCount;
	g_MuscleWeights = new WEIGHT[uLeafCount];
	CalcClustalWWeights(tree, g_MuscleWeights);
	}

void SetClustalWWeightsMuscle(MSA &msa)
	{
	if (0 == g_MuscleWeights)
		Quit("g_MuscleWeights = 0");
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const unsigned uId = msa.GetSeqId(uSeqIndex);
		if (uId >= g_uMuscleIdCount)
			Quit("SetClustalWWeightsMuscle: id out of range");
		msa.SetSeqWeight(uSeqIndex, g_MuscleWeights[uId]);
		}
	msa.NormalizeWeights((WEIGHT) 1.0);
	}

#define	LOCAL_VERBOSE	0

void SetThreeWayWeightsMuscle(MSA &msa)
	{
	if (NULL_NEIGHBOR == g_uTreeSplitNode1 || NULL_NEIGHBOR == g_uTreeSplitNode2)
		{
		msa.SetHenikoffWeightsPB();
		return;
		}

	const unsigned uMuscleSeqCount = g_ptrMuscleTree->GetLeafCount();
	WEIGHT *Weights = new WEIGHT[uMuscleSeqCount];

	CalcThreeWayWeights(*g_ptrMuscleTree, g_uTreeSplitNode1, g_uTreeSplitNode2,
	  Weights);

	const unsigned uMSASeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uMSASeqCount; ++uSeqIndex)
		{
		const unsigned uId = msa.GetSeqId(uSeqIndex);
		if (uId >= uMuscleSeqCount)
			Quit("SetThreeWayWeightsMuscle: id out of range");
		msa.SetSeqWeight(uSeqIndex, Weights[uId]);
		}
#if	LOCAL_VERBOSE
	{
	Log("SetThreeWayWeightsMuscle\n");
	for (unsigned n = 0; n < uMSASeqCount; ++n)
		{
		const unsigned uId = msa.GetSeqId(n);
		Log("%20.20s %6.3f\n", msa.GetSeqName(n), Weights[uId]);
		}
	}
#endif
	msa.NormalizeWeights((WEIGHT) 1.0);

	delete[] Weights;
	}

// Append msa2 at the end of msa1
void MSAAppend(MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount = msa1.GetSeqCount();

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	const unsigned uColCountCat = uColCount1 + uColCount2;

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
			{
			const char c = msa2.GetChar(uSeqIndex2, uColIndex);
			msa1.SetChar(uSeqIndex, uColCount1 + uColIndex, c);
			}
		}
	}

// "Catenate" two MSAs (by bad analogy with UNIX cat command).
// msa1 and msa2 must have same sequence names, but possibly
// in a different order.
// msaCat is the combined alignment produce by appending
// sequences in msa2 to sequences in msa1.
void MSACat(const MSA &msa1, const MSA &msa2, MSA &msaCat)
	{
	const unsigned uSeqCount = msa1.GetSeqCount();

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	const unsigned uColCountCat = uColCount1 + uColCount2;

	msaCat.SetSize(uSeqCount, uColCountCat);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		for (unsigned uColIndex = 0; uColIndex < uColCount1; ++uColIndex)
			{
			const char c = msa1.GetChar(uSeqIndex, uColIndex);
			msaCat.SetChar(uSeqIndex, uColIndex, c);
			}

		const char *ptrSeqName = msa1.GetSeqName(uSeqIndex);
		unsigned uSeqIndex2;
		msaCat.SetSeqName(uSeqIndex, ptrSeqName);
		bool bFound = msa2.GetSeqIndex(ptrSeqName, &uSeqIndex2);
		if (bFound)
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				{
				const char c = msa2.GetChar(uSeqIndex2, uColIndex);
				msaCat.SetChar(uSeqIndex, uColCount1 + uColIndex, c);
				}
			}
		else
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				msaCat.SetChar(uSeqIndex, uColCount1 + uColIndex, '-');
			}
		}
	}
