#include "muscle.h"
#include "profile.h"
#include "diaglist.h"

#define TRACE	0

#define pow4(i)	(1 << (2*i))	// 4^i = 2^(2*i)
const unsigned K = 7;
const unsigned KTUPS = pow4(K);
static unsigned TuplePos[KTUPS];

static char *TupleToStr(int t)
	{
	static char s[K];

	for (int i = 0; i < K; ++i)
		{
		unsigned Letter = (t/(pow4(i)))%4;
		assert(Letter >= 0 && Letter < 4);
		s[K-i-1] = LetterToChar(Letter);
		}

	return s;
	}

static unsigned GetTuple(const ProfPos *PP, unsigned uPos)
	{
	unsigned t = 0;

	for (unsigned i = 0; i < K; ++i)
		{
		const unsigned uLetter = PP[uPos+i].m_uResidueGroup;
		if (RESIDUE_GROUP_MULTIPLE == uLetter)
			return EMPTY;
		t = t*4 + uLetter;
		}

	return t;
	}

void FindDiagsNuc(const ProfPos *PX, unsigned uLengthX, const ProfPos *PY,
  unsigned uLengthY, DiagList &DL)
	{
	if (ALPHA_DNA != g_Alpha && ALPHA_RNA != g_Alpha)
		Quit("FindDiagsNuc: requires nucleo alphabet");

	DL.Clear();

// 16 is arbitrary slop, no principled reason for this.
	if (uLengthX < K + 16 || uLengthY < K + 16)
		return;

// Set A to shorter profile, B to longer
	const ProfPos *PA;
	const ProfPos *PB;
	unsigned uLengthA;
	unsigned uLengthB;
	bool bSwap;
	if (uLengthX < uLengthY)
		{
		bSwap = false;
		PA = PX;
		PB = PY;
		uLengthA = uLengthX;
		uLengthB = uLengthY;
		}
	else
		{
		bSwap = true;
		PA = PY;
		PB = PX;
		uLengthA = uLengthY;
		uLengthB = uLengthX;
		}

#if	TRACE
	Log("FindDiagsNuc(LengthA=%d LengthB=%d\n", uLengthA, uLengthB);
#endif

// Build tuple map for the longer profile, B
	if (uLengthB < K)
		Quit("FindDiags: profile too int");

	memset(TuplePos, EMPTY, sizeof(TuplePos));

	for (unsigned uPos = 0; uPos < uLengthB - K; ++uPos)
		{
		const unsigned uTuple = GetTuple(PB, uPos);
		if (EMPTY == uTuple)
			continue;
		TuplePos[uTuple] = uPos;
		}

// Find matches
	for (unsigned uPosA = 0; uPosA < uLengthA - K; ++uPosA)
		{
		const unsigned uTuple = GetTuple(PA, uPosA);
		if (EMPTY == uTuple)
			continue;
		const unsigned uPosB = TuplePos[uTuple];
		if (EMPTY == uPosB)
			continue;

	// This tuple is found in both profiles
		unsigned uStartPosA = uPosA;
		unsigned uStartPosB = uPosB;

	// Try to extend the match forwards
		unsigned uEndPosA = uPosA + K - 1;
		unsigned uEndPosB = uPosB + K - 1;
		for (;;)
			{
			if (uLengthA - 1 == uEndPosA || uLengthB - 1 == uEndPosB)
				break;
			const unsigned uAAGroupA = PA[uEndPosA+1].m_uResidueGroup;
			if (RESIDUE_GROUP_MULTIPLE == uAAGroupA)
				break;
			const unsigned uAAGroupB = PB[uEndPosB+1].m_uResidueGroup;
			if (RESIDUE_GROUP_MULTIPLE == uAAGroupB)
				break;
			if (uAAGroupA != uAAGroupB)
				break;
			++uEndPosA;
			++uEndPosB;
			}
		uPosA = uEndPosA;

#if	TRACE
		{
		Log("Match: A %4u-%4u   ", uStartPosA, uEndPosA);
		for (unsigned n = uStartPosA; n <= uEndPosA; ++n)
			Log("%c", LetterToChar(PA[n].m_uResidueGroup));
		Log("\n");
		Log("       B %4u-%4u   ", uStartPosB, uEndPosB);
		for (unsigned n = uStartPosB; n <= uEndPosB; ++n)
			Log("%c", LetterToChar(PB[n].m_uResidueGroup));
		Log("\n");
		}
#endif

		const unsigned uLength = uEndPosA - uStartPosA + 1;
		assert(uEndPosB - uStartPosB + 1 == uLength);

		if (uLength >= g_uMinDiagLength)
			{
			if (bSwap)
				DL.Add(uStartPosB, uStartPosA, uLength);
			else
				DL.Add(uStartPosA, uStartPosB, uLength);
			}
		}
	}
