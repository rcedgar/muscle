#include "muscle.h"
#include "mega.h"

void cmd_test_mega()
	{
#if 0
	Mega M;
	M.FromFile(g_Arg1);
	asserta(SIZE(M.m_Seqs) >= 2);
	uint index_X = 0;
	uint index_Y = 1;

	for(uint i = 0; i < (uint)M.m_Labels.size(); i++)
		{
		if(M.m_Labels[i] == "1hhs_A")
			index_X = i;
		if(M.m_Labels[i] == "1ra6_A")
			index_Y = i;
		//        cout << M.m_Labels[i] << "\t" << M.m_Seqs[i].size() << endl;
		}
	SetAlphaLC(false);

	string PWPath;
	float ea = AlignPairFlat_mega(&M, PWPath, index_X, index_Y);


	Sequence *InputSeq = Sequence::_NewSequence();
	InputSeq->FromString(M.m_Labels[index_X], M.m_Seqs[index_X]);
	Sequence *RefSeq = Sequence::_NewSequence();
	RefSeq->FromString(M.m_Labels[index_Y], M.m_Seqs[index_Y]);

	LogAln(*InputSeq, *RefSeq, PWPath);

	Sequence::_DeleteSequence(InputSeq);
	Sequence::_DeleteSequence(RefSeq);
#endif
	}
