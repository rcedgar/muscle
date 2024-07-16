#include "muscle.h"

void LoadInput(MultiSequence &InputSeqs)
	{
	if (opt(mega) || EndsWith(g_Arg1, ".mega"))
		{
		Mega::FromFile(g_Arg1);
		InputSeqs.FromStrings(Mega::m_Labels, Mega::m_Seqs);
		}
	else
		InputSeqs.LoadMFA(g_Arg1, true);
	SetGlobalInputMS(InputSeqs);
	}
