#include "muscle.h"
#include "profile.h"

void SetTermGaps(const ProfPos *Prof, unsigned uLength)
	{
	if (0 == uLength)
		return;

	ProfPos *First = (ProfPos *) Prof;
	ProfPos *Last = (ProfPos *) (Prof + uLength - 1);

	switch (g_TermGaps)
		{
	case TERMGAPS_Full:
		break;

	case TERMGAPS_Half:
	// -infinity check for lock left/right
		if (First->m_scoreGapOpen != MINUS_INFINITY)
			First->m_scoreGapOpen = 0;

		if (uLength > 1 && Last->m_scoreGapClose != MINUS_INFINITY)
			Last->m_scoreGapClose = 0;

	case TERMGAPS_Ext:
		if (First->m_scoreGapOpen != MINUS_INFINITY)
			First->m_scoreGapOpen *= -1;

		if (uLength > 1 && Last->m_scoreGapClose != MINUS_INFINITY)
			Last->m_scoreGapClose *= -1;
		break;

	default:
		Quit("Invalid g_TermGaps");
		}
	}
