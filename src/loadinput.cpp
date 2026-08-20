#include "muscle.h"
#include "mega.h"

void LoadInput(MultiSequence &InputSeqs)
	{
	const bool WantMega = opt(mega) || EndsWith(g_Arg1, ".mega");
	const bool WantStructs = opt(structs) || Mega::IsStructsInput(g_Arg1);

	if (WantMega && WantStructs)
		Die("Input looks like both .mega and STRUCTS; use one");

	if (WantMega)
		{
		Mega::FromFile(g_Arg1);
		InputSeqs.FromStrings(Mega::m_Labels, Mega::m_Seqs);
		}
	else if (WantStructs)
		{
		Mega::FromStructs(g_Arg1);
		InputSeqs.FromStrings(Mega::m_Labels, Mega::m_Seqs);
		}
	else
		InputSeqs.LoadMFA(g_Arg1, true);
	SetGlobalInputMS(InputSeqs);
	}
