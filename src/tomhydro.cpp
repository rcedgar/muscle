#include "muscle.h"
#include "profile.h"

// Original:
//HYDROPHILIC_CONTEXT 0 6 -0.3969495574
//HYDROPHILIC_CONTEXT 1 6 -0.9407126603
//HYDROPHILIC_CONTEXT 2 6 -0.4968150972
//HYDROPHILIC_CONTEXT 3 6 -0.271646023
//HYDROPHILIC_CONTEXT 4 6 0.006990406416
//HYDROPHILIC_CONTEXT 5 6 0.1381111256
//HYDROPHILIC_CONTEXT 6 6 0.2541439872

// Blosum62:
//HYDROPHILIC_CONTEXT 0 6 -0.2448419585
//HYDROPHILIC_CONTEXT 1 6 -0.8734889946
//HYDROPHILIC_CONTEXT 2 6 -0.5724336598
//HYDROPHILIC_CONTEXT 3 6 -0.2670439975
//HYDROPHILIC_CONTEXT 4 6 0.004844647323
//HYDROPHILIC_CONTEXT 5 6 0.1812057148
//HYDROPHILIC_CONTEXT 6 6 0.1036540864

static SCORE Factors[7] =
	{
	(SCORE) -0.2448419585,
	(SCORE) -0.8734889946,
	(SCORE) -0.5724336598,
	(SCORE) -0.2670439975,
	(SCORE) 0.004844647323,
	(SCORE) 0.1812057148,
	(SCORE) 0.1036540864
	};

static bool Hydrophilic[20] =
	{
	false,		// A
	false,		// C
	true,		// D
	true,		// E
	false,		// F
	true,		// G
	false,		// H
	false,		// I
	true,		// K
	false,		// L
	false,		// M
	true,		// N
	true,		// P
	true,		// Q
	true,		// R
	true,		// S
	false,		// T
	false,		// V
	false,		// Y
	false,		// W
	};

bool IsHydrophilic(const FCOUNT fcCounts[])
	{
	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		if (fcCounts[uLetter] > 0.0 && Hydrophilic[uLetter])
			return false;
	return true;
	}

static double HydrophilicFraction(const FCOUNT fcCounts[])
	{
	double TotalAll = 0.0;
	double TotalHydrophilic = 0.0;
	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		{
		FCOUNT Freq = fcCounts[uLetter];
		TotalAll += Freq;
		if (Hydrophilic[uLetter])
			TotalHydrophilic += Freq;
		}
	return TotalHydrophilic / TotalAll;
	}

void TomHydro(ProfPos *Prof, unsigned uLength)
	{
	if (ALPHA_Amino != g_Alpha)
		return;
	if (uLength < 6)
		return;

	for (unsigned uColIndex = 3; uColIndex < uLength - 2; ++uColIndex)
		{
	// 6-residue window:
	//	 xxxxxx
	//	AARNCARNGTAGCATNAC
	//	AARN----------TNAC

		double dCount = 0.0;
		for (unsigned uColIndexW = uColIndex - 3; uColIndexW < uColIndex + 3;
		  ++uColIndexW)
			{
			const ProfPos &PP = Prof[uColIndexW];
			dCount += HydrophilicFraction(PP.m_fcCounts);
			}
	// Round to nearest integer
		unsigned uCount = (unsigned) (dCount + 0.5);
		if (uCount > 6)
			uCount = 6;
		SCORE dFactor = Factors[uCount];
		ProfPos &PP = Prof[uColIndex];
		PP.m_scoreGapOpen += dFactor;
		PP.m_scoreGapClose += dFactor;
		}
	}
