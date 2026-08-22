#include "muscle.h"
#include "mega.h"

void cmd_test_mega()
	{
#if 0
	Mega::FromStructs(g_Arg1);
	asserta(SIZE(Mega::m_Seqs) >= 2);
	uint index_X = 0;
	uint index_Y = 1;

	for(uint i = 0; i < (uint)Mega::m_Labels.size(); i++)
		{
		if(Mega::m_Labels[i] == "1hhs_A")
			index_X = i;
		if(Mega::m_Labels[i] == "1ra6_A")
			index_Y = i;
		}
	SetAlphaLC(false);

	string PWPath;
	float ea = AlignPairFlat_mega(0, PWPath, index_X, index_Y);

	Sequence *InputSeq = Sequence::_NewSequence();
	InputSeq->FromString(Mega::m_Labels[index_X], Mega::m_Seqs[index_X]);
	Sequence *RefSeq = Sequence::_NewSequence();
	RefSeq->FromString(Mega::m_Labels[index_Y], Mega::m_Seqs[index_Y]);

	LogAln(*InputSeq, *RefSeq, PWPath);

	Sequence::_DeleteSequence(InputSeq);
	Sequence::_DeleteSequence(RefSeq);
#endif
	}
