#include "muscle.h"

//// Pascaralle and Argos gap factors
//// after Table 1 in Thompson et. al. ClustalW NAR paper.
//static double PAFFacs[20] =
//	{
//	1.13,		// A
//	1.13,		// C
//	0.96,		// D
//	1.31,		// E
//	1.20,		// F
//	0.61,		// G
//	1.00,		// H
//	1.32,		// I
//	0.96,		// K
//	1.21,		// L
//	1.29,		// M
//	0.62,		// N
//	0.74,		// P
//	1.07,		// Q
//	0.72,		// R
//	0.76,		// S
//	0.89,		// T
//	1.25,		// V
//	1.00,		// Y
//	1.23,		// W
//	};
//
//// (Not used: does not appear to work well).
//SCORE PAFactor(const FCOUNT fcCounts[])
//	{
//	if (ALPHA_Amino != g_Alpha)
//		Quit("PAFFactor: requires amino acid sequence");
//
//	FCOUNT fLetterCount = 0;
//	double dSum = 0;
//	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
//		{
//		const FCOUNT fCount = fcCounts[uLetter];
//		dSum += fCount*PAFFacs[uLetter];
//		fLetterCount += fCount;
//		}
//	if (0 == fLetterCount)
//		return 0.5;
//	return (SCORE) (dSum/fLetterCount);
//	}

//static bool Hydrophilic[20] =
//	{
//	false,		// A
//	false,		// C
//	true,		// D
//	true,		// E
//	false,		// F
//	true,		// G
//	false,		// H
//	false,		// I
//	true,		// K
//	false,		// L
//	false,		// M
//	true,		// N
//	true,		// P
//	true,		// Q
//	true,		// R
//	true,		// S
//	false,		// T
//	false,		// V
//	false,		// Y
//	false,		// W
//	};
//
//bool IsHydrophilic(const FCOUNT fcCounts[])
//	{
//	if (ALPHA_Amino != g_Alpha)
//		Quit("IsHydrophilic: requires amino acid sequence");
//
//	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
//		if (fcCounts[uLetter] > 0 && !Hydrophilic[uLetter])
//			return false;
//	return true;
//	}
//
//bool IsHydrophilic(const unsigned uCounts[])
//	{
//	if (ALPHA_Amino != g_Alpha)
//		Quit("IsHydrophilic: requires amino acid sequence");
//
//	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
//		if (uCounts[uLetter] > 0 && !Hydrophilic[uLetter])
//			return false;
//	return true;
//	}

// LIVCATMFYWHK
// Venn		Pascaralla	B&T		Me
// L		y			y		y
// I		y			y		y
// V		y			y		y
// C		y			n
// A		y			y		y
// T		N			n
// M		y			y		y
// F		y			y		y
// Y		n			n
// W		y			n
// H		n			n
// K		n			n
static bool Hydrophobic[20] =
	{
	true,		// A
	true,		// C
	false,		// D
	false,		// E
	true,		// F
	false,		// G
	true,		// H
	true,		// I
	false,		// K
	true,		// L
	true,		// M
	false,		// N
	false,		// P
	false,		// Q
	false,		// R
	false,		// S
	true,		// T
	true,		// V
	true,		// Y
	true,		// W
	};

bool IsHydrophobic(const FCOUNT fcCounts[])
	{
	if (ALPHA_Amino != g_Alpha)
		Quit("IsHydrophobic: requires amino acid sequence");

	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		if (fcCounts[uLetter] > 0.0 && !Hydrophobic[uLetter])
			return false;
	return true;
	}
