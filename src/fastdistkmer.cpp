#include "muscle.h"
#include "msa.h"
#include "seqvect.h"
#include "seq.h"
#include "distfunc.h"
#include <math.h>

#define TRACE	0

/***
Some candidate alphabets considered because they
have high correlations and small table sizes.
Correlation coefficent is between k-mer distance
and %id D measured from a CLUSTALW alignment.
Table size is N^k where N is size of alphabet.
A is standard (uncompressed) amino alphabet.

                           Correlation
Alpha   N  k  Table Size   all   25-50%
-----  --  -  ----------   ----  ------
A      20  3       8,000  0.943   0.575
A      20  4     160,000  0.962   0.685 <<
LiA    14  4      38,416  0.966   0.645
SEB    14  4      38,416  0.964   0.634
LiA    13  4      28,561  0.965   0.640
LiA    12  4      20,736  0.963   0.620
LiA    10  5     100,000  0.964   0.652

We select A with k=4 because it has the best
correlations. The only drawback is a large table
size, but space is readily available and the only 
additional time cost is in resetting the table to
zero, which can be done quickly with memset or by
keeping a list of the k-mers that were found (should
test to see which is faster, and may vary by compiler
and processor type). It also has the minor advantage
that we don't need to convert the alphabet.

Fractional identity d is estimated as follows.

	F = fractional k-mer count
	if F is 0: F = 0.01
	Y = log(0.02 + F)
	d = -4.1 + 4.12*Y

The constant 0.02 was chosen to make the relationship
between Y and D linear. The constants -4.1 and 4.12
were chosen to fit a straight line to the scatterplot
of Y vs D.
***/

#define MIN(x, y)	(((x) < (y)) ? (x) : (y))

const unsigned K = 4;
const unsigned N = 20;
const unsigned N_2 = 20*20;
const unsigned N_3 = 20*20*20;
const unsigned N_4 = 20*20*20*20;

const unsigned TABLE_SIZE = N_4;

// For debug output
const char *KmerToStr(unsigned Kmer)
	{
	static char s[5];

	unsigned c3 = (Kmer/N_3)%N;
	unsigned c2 = (Kmer/N_2)%N;
	unsigned c1 = (Kmer/N)%N;
	unsigned c0 = Kmer%N;

	s[0] = LetterToChar(c3);
	s[1] = LetterToChar(c2);
	s[2] = LetterToChar(c1);
	s[3] = LetterToChar(c0);
	return s;
	}

void CountKmers(const byte s[], unsigned uSeqLength, byte KmerCounts[])
	{
#if	TRACE
	Log("CountKmers\n");
#endif
	memset(KmerCounts, 0, TABLE_SIZE*sizeof(byte));

	const byte *ptrKmerStart = s;
	const byte *ptrKmerEnd = s + 4;
	const byte *ptrSeqEnd = s + uSeqLength;

	unsigned c3 = s[0]*N_3;
	unsigned c2 = s[1]*N_2;
	unsigned c1 = s[2]*N;
	unsigned c0 = s[3];

	unsigned Kmer = c3 + c2 + c1 + c0;

	for (;;)
		{
		assert(Kmer < TABLE_SIZE);

#if	TRACE
		Log("Kmer=%d=%s\n", Kmer, KmerToStr(Kmer));
#endif
		++(KmerCounts[Kmer]);

		if (ptrKmerEnd == ptrSeqEnd)
			break;

	// Compute k-mer as function of previous k-mer:
	// 1. Subtract first letter from previous k-mer.
	// 2. Multiply by N.
	// 3. Add next letter.
		c3 = (*ptrKmerStart++) * N_3;
		Kmer = (Kmer - c3)*N;
		Kmer += *ptrKmerEnd++;
		}
	}

unsigned CommonKmerCount(const byte Seq[], unsigned uSeqLength,
  const byte KmerCounts1[], const byte Seq2[], unsigned uSeqLength2)
	{
	byte KmerCounts2[TABLE_SIZE];
	CountKmers(Seq2, uSeqLength2, KmerCounts2);

	const byte *ptrKmerStart = Seq;
	const byte *ptrKmerEnd = Seq + 4;
	const byte *ptrSeqEnd = Seq + uSeqLength;

	unsigned c3 = Seq[0]*N_3;
	unsigned c2 = Seq[1]*N_2;
	unsigned c1 = Seq[2]*N;
	unsigned c0 = Seq[3];

	unsigned Kmer = c3 + c2 + c1 + c0;

	unsigned uCommonCount = 0;
	for (;;)
		{
		assert(Kmer < TABLE_SIZE);

		const byte Count1 = KmerCounts1[Kmer];
		const byte Count2 = KmerCounts2[Kmer];

		uCommonCount += MIN(Count1, Count2);

	// Hack so we don't double-count
		KmerCounts2[Kmer] = 0;

		if (ptrKmerEnd == ptrSeqEnd)
			break;

	// Compute k-mer as function of previous k-mer:
	// 1. Subtract first letter from previous k-mer.
	// 2. Multiply by N.
	// 3. Add next letter.
		c3 = (*ptrKmerStart++) * N_3;
		Kmer = (Kmer - c3)*N;
		Kmer += *ptrKmerEnd++;
		}
	return uCommonCount;
	}

static void SeqToLetters(const Seq &s, byte Letters[])
	{
	const unsigned uSeqLength = s.Length();
	for (unsigned uCol = 0; uCol < uSeqLength; ++uCol)
		{
		char c = s.GetChar(uCol);
	// Ugly hack. My k-mer counting code isn't wild-card
	// aware. Arbitrarily replace wildcards by a specific
	// amino acid.
		if (IsWildcardChar(c))
			c = 'A';
		*Letters++ = CharToLetter(c);
		}
	}

void FastDistKmer(const SeqVect &v, DistFunc &DF)
	{
	byte KmerCounts[TABLE_SIZE];

	const unsigned uSeqCount = v.GetSeqCount();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;

// Initialize distance matrix to zero
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

	unsigned uMaxLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq &s = v.GetSeq(uSeqIndex);
		unsigned uSeqLength = s.Length();
		if (uSeqLength > uMaxLength)
			uMaxLength = uSeqLength;
		}
	if (0 == uMaxLength)
		return;

	byte *Seq1Letters = new byte[uMaxLength];
	byte *Seq2Letters = new byte[uMaxLength];

	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount - 1; ++uSeqIndex1)
		{
		const Seq &s1 = v.GetSeq(uSeqIndex1);
		const unsigned uSeqLength1 = s1.Length();

		SeqToLetters(s1, Seq1Letters);
		CountKmers(Seq1Letters, uSeqLength1, KmerCounts);

		for (unsigned uSeqIndex2 = uSeqIndex1 + 1; uSeqIndex2 < uSeqCount;
		  ++uSeqIndex2)
			{
			const Seq &s2 = v.GetSeq(uSeqIndex2);
			const unsigned uSeqLength2 = s2.Length();

			SeqToLetters(s2, Seq2Letters);

			unsigned uCommonKmerCount = CommonKmerCount(Seq1Letters, uSeqLength1,
			  KmerCounts, Seq2Letters, uSeqLength2);

			unsigned uMinLength = MIN(uSeqLength1, uSeqLength2);
			double F = (double) uCommonKmerCount / (uMinLength - K + 1);
			if (0.0 == F)
				F = 0.01;
			double Y = log(0.02 + F);
			double EstimatedPctId = Y/4.12 + 0.995;
			double KD = KimuraDist(EstimatedPctId);
//			DF.SetDist(uSeqIndex1, uSeqIndex2, (float) KD);
			DF.SetDist(uSeqIndex1, uSeqIndex2, (float) (1 - F));
#if	TRACE
			Log("CommonCount=%u, MinLength=%u, F=%6.4f Y=%6.4f, %%id=%6.4f, KimuraDist=%8.4f\n",
			  uCommonKmerCount, uMinLength, F, Y, EstimatedPctId, KD);
#endif
			}
		}

	delete[] Seq1Letters;
	delete[] Seq2Letters;
	}
