#include "muscle.h"
#include "msa.h"

static int Blosum62[23][23] =
	{
//   A   B   C   D   E    F   G   H   I   K    L   M   N   P   Q    R   S   T   V   W    X   Y   Z 
	+4, -2, +0, -2, -1,  -2, +0, -2, -1, -1,  -1, -1, -2, -1, -1,  -1, +1, +0, +0, -3,  -1, -2, -1,  // A
	-2, +6, -3, +6, +2,  -3, -1, -1, -3, -1,  -4, -3, +1, -1, +0,  -2, +0, -1, -3, -4,  -1, -3, +2,  // B
	+0, -3, +9, -3, -4,  -2, -3, -3, -1, -3,  -1, -1, -3, -3, -3,  -3, -1, -1, -1, -2,  -1, -2, -4,  // C
	-2, +6, -3, +6, +2,  -3, -1, -1, -3, -1,  -4, -3, +1, -1, +0,  -2, +0, -1, -3, -4,  -1, -3, +2,  // D
	-1, +2, -4, +2, +5,  -3, -2, +0, -3, +1,  -3, -2, +0, -1, +2,  +0, +0, -1, -2, -3,  -1, -2, +5,  // E
	
	-2, -3, -2, -3, -3,  +6, -3, -1, +0, -3,  +0, +0, -3, -4, -3,  -3, -2, -2, -1, +1,  -1, +3, -3,  // F
	+0, -1, -3, -1, -2,  -3, +6, -2, -4, -2,  -4, -3, +0, -2, -2,  -2, +0, -2, -3, -2,  -1, -3, -2,  // G
	-2, -1, -3, -1, +0,  -1, -2, +8, -3, -1,  -3, -2, +1, -2, +0,  +0, -1, -2, -3, -2,  -1, +2, +0,  // H
	-1, -3, -1, -3, -3,  +0, -4, -3, +4, -3,  +2, +1, -3, -3, -3,  -3, -2, -1, +3, -3,  -1, -1, -3,  // I
	-1, -1, -3, -1, +1,  -3, -2, -1, -3, +5,  -2, -1, +0, -1, +1,  +2, +0, -1, -2, -3,  -1, -2, +1,  // K
	
	-1, -4, -1, -4, -3,  +0, -4, -3, +2, -2,  +4, +2, -3, -3, -2,  -2, -2, -1, +1, -2,  -1, -1, -3,  // L
	-1, -3, -1, -3, -2,  +0, -3, -2, +1, -1,  +2, +5, -2, -2, +0,  -1, -1, -1, +1, -1,  -1, -1, -2,  // M
	-2, +1, -3, +1, +0,  -3, +0, +1, -3, +0,  -3, -2, +6, -2, +0,  +0, +1, +0, -3, -4,  -1, -2, +0,  // N
	-1, -1, -3, -1, -1,  -4, -2, -2, -3, -1,  -3, -2, -2, +7, -1,  -2, -1, -1, -2, -4,  -1, -3, -1,  // P
	-1, +0, -3, +0, +2,  -3, -2, +0, -3, +1,  -2, +0, +0, -1, +5,  +1, +0, -1, -2, -2,  -1, -1, +2,  // Q
	
	-1, -2, -3, -2, +0,  -3, -2, +0, -3, +2,  -2, -1, +0, -2, +1,  +5, -1, -1, -3, -3,  -1, -2, +0,  // R
	+1, +0, -1, +0, +0,  -2, +0, -1, -2, +0,  -2, -1, +1, -1, +0,  -1, +4, +1, -2, -3,  -1, -2, +0,  // S
	+0, -1, -1, -1, -1,  -2, -2, -2, -1, -1,  -1, -1, +0, -1, -1,  -1, +1, +5, +0, -2,  -1, -2, -1,  // T
	+0, -3, -1, -3, -2,  -1, -3, -3, +3, -2,  +1, +1, -3, -2, -2,  -3, -2, +0, +4, -3,  -1, -1, -2,  // V
	-3, -4, -2, -4, -3,  +1, -2, -2, -3, -3,  -2, -1, -4, -4, -2,  -3, -3, -2, -3,+11,  -1, +2, -3,  // W
	
	-1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1, -1, -1,  -1, -1, -1,  // X
	-2, -3, -2, -3, -2,  +3, -3, +2, -1, -2,  -1, -1, -2, -3, -1,  -2, -2, -2, -1, +2,  -1, +7, -2,  // Y
	-1, +2, -4, +2, +5,  -3, -2, +0, -3, +1,  -3, -2, +0, -1, +2,  +0, +0, -1, -2, -3,  -1, -2, +5,  // Z
	};

static int toi_tab[26] =
	{
	0,	// A
	1,	// B
	2,	// C
	3,	// D
	4,	// E
	5,	// F
	6,	// G
	7,	// H
	8,	// I
	-1,	// J
	9,	// K
	10,	// L
	11,	// M
	12,	// N
	-1,	// O
	13,	// P
	14,	// Q
	15,	// R
	16,	// S
	17,	// T
	17,	// U
	18,	// V
	19,	// W
	20,	// X
	21,	// Y
	22,	// Z
	};

static int toi(char c)
	{
	c = toupper(c);
	return toi_tab[c - 'A'];
	}

static int BlosumScore(char c1, char c2)
	{
	int i1 = toi(c1);
	int i2 = toi(c2);
	return Blosum62[i1][i2];
	}

/***
Consider a column with 5 As and 3 Bs.
There are:
	5x4 pairs of As.
	3x2 pairs of Bs.
	5x3x2 AB pairs
	8x7 = 5x4 + 3x2 + 5x3x2 pairs of letters
***/
static double BlosumScoreCol(const MSA &a, unsigned uColIndex)
	{
	int iCounts[23];
	memset(iCounts, 0, sizeof(iCounts));
	const unsigned uSeqCount = a.GetSeqCount();
	unsigned uCharCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		char c = a.GetChar(uSeqIndex, uColIndex);
		if (IsGapChar(c))
			continue;
		int iChar = toi(c);
		++iCounts[iChar];
		++uCharCount;
		}
	if (uCharCount < 2)
		return -9;
	int iTotalScore = 0;
	for (int i1 = 0; i1 < 23; ++i1)
		{
		int iCounts1 = iCounts[i1];
		iTotalScore += iCounts1*(iCounts1 - 1)*Blosum62[i1][i1];
		for (int i2 = i1 + 1; i2 < 23; ++i2)
			iTotalScore += iCounts[i2]*iCounts1*2*Blosum62[i1][i2];
		}
	int iPairCount = uCharCount*(uCharCount - 1);
	return (double) iTotalScore / (double) iPairCount;
	}

/***
Consider a column with 5 As and 3 Bs.
A residue of type Q scores:
	5xAQ + 3xBQ
***/
static void AssignColorsCol(const MSA &a, unsigned uColIndex, int **Colors)
	{
	int iCounts[23];
	memset(iCounts, 0, sizeof(iCounts));
	const unsigned uSeqCount = a.GetSeqCount();
	unsigned uCharCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		char c = a.GetChar(uSeqIndex, uColIndex);
		if (IsGapChar(c))
			continue;
		int iChar = toi(c);
		++iCounts[iChar];
		++uCharCount;
		}
	int iMostConservedType = -1;
	int iMostConservedCount = -1;
	for (unsigned i = 0; i < 23; ++i)
		{
		if (iCounts[i] > iMostConservedCount)
			{
			iMostConservedType = i;
			iMostConservedCount = iCounts[i];
			}
		}

	double dColScore = BlosumScoreCol(a, uColIndex);
	int c;
	if (dColScore >= 3.0)
		c = 3;
	//else if (dColScore >= 1.0)
	//	c = 2;
	else if (dColScore >= 0.2)
		c = 1;
	else
		c = 0;

	int Color[23];
	for (unsigned uLetter = 0; uLetter < 23; ++uLetter)
		{
		double dScore = Blosum62[uLetter][iMostConservedType];
		if (dScore >= dColScore)
			Color[uLetter] = c;
		else
			Color[uLetter] = 0;
		}

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		char c = a.GetChar(uSeqIndex, uColIndex);
		if (IsGapChar(c))
			{
			Colors[uSeqIndex][uColIndex] = 0;
			continue;
			}
		int iLetter = toi(c);
		if (iLetter >= 0 && iLetter < 23)
			Colors[uSeqIndex][uColIndex] = Color[iLetter];
		else
			Colors[uSeqIndex][uColIndex] = 0;
		}
	}

void AssignColors(const MSA &a, int **Colors)
	{
	const unsigned uColCount = a.GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		AssignColorsCol(a, uColIndex, Colors);
	}
