#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"
#include <math.h>

#define	MIN(x, y)	((x) < (y) ? (x) : (y))

static void SetKmerBitVector(const Seq &s, byte Bits[])
	{
	const unsigned uLength = s.Length();
	const unsigned k = 3;	// kmer length
	unsigned i = 0;
	unsigned c = 0;
	unsigned h = 0;
	for (unsigned j = 0; j < k - 1; ++j)
		{
		unsigned x = CharToLetterEx(s[i++]);
		if (x <= AX_Y)
			c = c*20 + x;
		else
			{
			c = 0;
			h = j + 1;
			}
		}
	for ( ; i < uLength; ++i)
		{
		unsigned x = CharToLetterEx(s[i++]);
		if (x <= AX_Y)
			c = (c*20 + x)%8000;
		else
			{
			c = 0;
			h = i + k;
			}
		if (i >= h)
			{
			unsigned ByteOffset = c/8;
			unsigned BitOffset = c%8;
			Bits[ByteOffset] |= (1 << BitOffset);
			}
		}
	}

static unsigned CommonBitCount(const byte Bits1[], const byte Bits2[])
	{
	const byte * const p1end = Bits1 + 1000;
	const byte *p2 = Bits2;

	unsigned uCount = 0;
	for (const byte *p1 = Bits1; p1 != p1end; ++p1)
		{
	// Here is a cute trick for efficiently counting the
	// bits common between two bytes by combining them into
	// a single word.
		unsigned b = *p1 | (*p2 << 8);
		while (b != 0)
			{
			if (b & 0x101)
				++uCount;
			b >>= 1;
			}
		++p2;
		}
	return uCount;
	}

void DistKbit20_3(const SeqVect &v, DistFunc &DF)
	{
	const unsigned uSeqCount = v.Length();
	DF.SetCount(uSeqCount);

// There are 20^3 = 8,000 distinct kmers in the 20-letter alphabet.
// For each sequence, we create a bit vector of length 8,000, i.e.
// 1,000 bytes, having one bit per kmer. The bit is set to 1 if the
// kmer is present in the sequence.
	const unsigned uBytes = uSeqCount*1000;
	byte *BitVector = new byte[uBytes];
	memset(BitVector, 0, uBytes);

	SetProgressDesc("K-bit distance matrix");
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		SetKmerBitVector(*v[uSeqIndex], BitVector + uSeqIndex*1000);

	unsigned uDone = 0;
	const unsigned uTotal = (uSeqCount*(uSeqCount - 1))/2;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const byte *Bits1 = BitVector + uSeqIndex1*1000;
		const unsigned uLength1 = v[uSeqIndex1]->Length();
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqIndex1; ++uSeqIndex2)
			{
			const byte *Bits2 = BitVector + uSeqIndex2*1000;
			const unsigned uLength2 = v[uSeqIndex2]->Length();
			const float fCount = (float) CommonBitCount(Bits1, Bits2);

		// Distance measure = K / min(L1, L2)
		// K is number of distinct kmers that are found in both sequences
			const float fDist = fCount / MIN(uLength1, uLength2);
			DF.SetDist(uSeqIndex1, uSeqIndex2, fDist);
			if (uDone%10000 == 0)
				Progress(uDone, uTotal);
			++uDone;
			}
		}
	ProgressStepsDone();

	delete[] BitVector;
	}
