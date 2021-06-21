#include "muscle.h"
#include "tree.h"
#include "profile.h"
#include "msa.h"
#include "seqvect.h"
#include "pwpath.h"

static void DoSeq(Seq &s, unsigned uSeqIndex, const ProfPos *RootProf,
  unsigned uRootProfLength, MSA &msaOut)
	{
	MSA msaSeq;
	msaSeq.FromSeq(s);
	const unsigned uSeqLength = s.Length();

	MSA msaDummy;
	msaDummy.SetSize(1, uRootProfLength);
	msaDummy.SetSeqId(0, 0);
	msaDummy.SetSeqName(0, "Dummy0");
	for (unsigned uColIndex = 0; uColIndex < uRootProfLength; ++uColIndex)
		msaDummy.SetChar(0, uColIndex, '?');

	ProfPos *SeqProf = ProfileFromMSA(msaSeq);
	for (unsigned uColIndex = 0; uColIndex < uSeqLength; ++uColIndex)
		{
		ProfPos &PP = SeqProf[uColIndex];
		PP.m_scoreGapOpen = MINUS_INFINITY;
		PP.m_scoreGapClose = MINUS_INFINITY;
		}

	ProfPos *ProfOut;
	unsigned uLengthOut;
	PWPath Path;
	AlignTwoProfs(SeqProf, uSeqLength, 1.0, RootProf, uRootProfLength, 1.0,
	  Path, &ProfOut, &uLengthOut);
	assert(uLengthOut = uRootProfLength);
	delete[] ProfOut;

	MSA msaCombined;
	AlignTwoMSAsGivenPath(Path, msaSeq, msaDummy, msaCombined);

	msaCombined.LogMe();
	msaOut.SetSeqName(uSeqIndex, s.GetName());
	msaOut.SetSeqId(uSeqIndex, s.GetId());
	for (unsigned uColIndex = 0; uColIndex < uRootProfLength; ++uColIndex)
		msaOut.SetChar(uSeqIndex, uColIndex, msaCombined.GetChar(0, uColIndex));
	}

// Steven Brenner's O(NL^2) proposal for creating a root alignment
// Align each sequence to the profile at the root.
// Compare the e-string solution, which is O(NL log N).
void MakeRootMSABrenner(SeqVect &v, const Tree &GuideTree, ProgNode Nodes[],
  MSA &a)
	{
	const unsigned uSeqCount = v.Length();
	const unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	const ProfPos *RootProfile = Nodes[uRootNodeIndex].m_Prof;
	const unsigned uRootColCount = Nodes[uRootNodeIndex].m_uLength;
	a.SetSize(uSeqCount, uRootColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		DoSeq(*v[uSeqIndex], uSeqIndex, RootProfile, uRootColCount, a);
	}
