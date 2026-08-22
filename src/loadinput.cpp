#include "muscle.h"
#include "mega.h"

void LoadInput(MultiSequence &InputSeqs)
	{
	Mega::RejectLegacyMega(g_Arg1);

	const bool WantStructs = opt(structs) || Mega::IsStructsInput(g_Arg1);

	if (WantStructs)
		{
		Mega::FromStructs(g_Arg1);
		InputSeqs.FromStrings(Mega::m_Labels, Mega::m_Seqs);
		}
	else
		InputSeqs.LoadMFA(g_Arg1, true);
	SetGlobalInputMS(InputSeqs);
	}
