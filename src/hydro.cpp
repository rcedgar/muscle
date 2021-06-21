#include "muscle.h"
#include "profile.h"

extern void TomHydro(ProfPos *Prof, unsigned Length);

// Apply hydrophobicity heuristic to a profile
void Hydro(ProfPos *Prof, unsigned uLength)
	{
	if (ALPHA_Amino != g_Alpha)
		return;

	if (g_bTomHydro)
		{
		TomHydro(Prof, uLength);
		return;
		}

	if (0 == g_uHydrophobicRunLength)
		return;

	if (uLength <= g_uHydrophobicRunLength)
		return;

	unsigned uRunLength = 0;
	unsigned L2 = g_uHydrophobicRunLength/2;
	for (unsigned uColIndex = L2; uColIndex < uLength - L2; ++uColIndex)
		{
		ProfPos &PP = Prof[uColIndex];
		bool bHydro = IsHydrophobic(PP.m_fcCounts);
		if (bHydro)
			{
			++uRunLength;
			if (uRunLength >= g_uHydrophobicRunLength)
				{
				Prof[uColIndex-L2].m_scoreGapOpen *= (SCORE) g_dHydroFactor;
				Prof[uColIndex-L2].m_scoreGapClose *= (SCORE) g_dHydroFactor;
				}
			}
		else
			uRunLength = 0;
		}
	}
