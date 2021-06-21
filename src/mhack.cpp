#include "muscle.h"
#include "seqvect.h"
#include "msa.h"

/***
Methionine hack.
Most proteins start with M.
This results in odd-looking alignments with the terminal Ms aligned followed
immediately by gaps.
Hack this by treating terminal M like X.
***/

static bool *M;

void MHackStart(SeqVect &v)
	{
	if (ALPHA_Amino != g_Alpha)
		return;

	const unsigned uSeqCount = v.Length();
	M = new bool[uSeqCount];
	memset(M, 0, uSeqCount*sizeof(bool));
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = v.GetSeq(uSeqIndex);
		if (0 == s.Length())
			continue;
		unsigned uId = s.GetId();
		if (s[0] == 'M' || s[0] == 'm')
			{
			M[uId] = true;
			s[0] = 'X';
			}
		}
	}

void MHackEnd(MSA &msa)
	{
	if (ALPHA_Amino != g_Alpha)
		return;
	if (0 == M)
		return;

	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uId = msa.GetSeqId(uSeqIndex);
		if (M[uId])
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
				{
				if (!msa.IsGap(uSeqIndex, uColIndex))
					{
					msa.SetChar(uSeqIndex, uColIndex, 'M');
					break;
					}
				}
			}
		}

	delete[] M;
	M = 0;
	}
