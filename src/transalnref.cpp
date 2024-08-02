#include "muscle.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "transaln.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path);

static void StripGaps(const string &Row, string &UngappedSeq)
	{
	for (uint i = 0; i < SIZE(Row); ++i)
		{
		byte r = Row[i];
		if (!isgap(r))
			UngappedSeq += r;
		}
	}

void cmd_transalnref()
	{
	asserta(optset_input2);
	asserta(optset_output);
	asserta(optset_label);

	const string &RefAlnFN = g_Arg1;
	const string &AddFastaFN = opt(input2);
	const string &RefLabel = opt(label);

	MultiSequence RefAln;
	RefAln.LoadMFA(RefAlnFN);
	const uint RefSeqCount = RefAln.GetSeqCount();
	uint RefSeqIndex = UINT_MAX;
	for (uint i = 0; i < RefSeqCount; ++i)
		{
		string Label = string(RefAln.GetLabel(i));
		if (Label == RefLabel)
			{
			RefSeqIndex = i;
			break;
			}
		}
	if (RefSeqIndex == UINT_MAX)
		Die("Not found >%s", RefLabel.c_str());

	const Sequence *RowR = RefAln.m_Seqs[RefSeqIndex];

	string RowRStr;
	RowR->GetSeqAsString(RowRStr);
	string RStr;
	StripGaps(RowRStr, RStr);
	const uint LR = SIZE(RStr);

	MultiSequence SeqsToAdd;
	SeqsToAdd.LoadMFA(AddFastaFN);
	if (SeqsToAdd.GetSeqCount() != 1)
		Die("-input2 must have exactly one sequence");
	const string AddLabel = string(SeqsToAdd.GetLabel(0));
	const Sequence *RowA = SeqsToAdd.m_Seqs[0];

	string RowAStr;
	RowA->GetSeqAsString(RowAStr);
	string AStr;
	StripGaps(RowAStr, AStr);
	const uint LA = SIZE(AStr);

	PathInfo *PI = ObjMgr::GetPathInfo();
	XDPMem Mem;
	ViterbiFastMem(Mem, (const byte *) RStr.c_str(), LR, (const byte *) AStr.c_str(), LA, *PI);
	string PathStr;
	PI->GetPathStr(PathStr);

	const uint PairColCount = PI->GetColCount();
	uint PosR = 0;
	uint PosA = 0;
	uint Ids = 0;
	for (uint Col = 0; Col < PairColCount; ++Col)
		{
		switch (PathStr[Col])
			{
		case 'M':
			{
			char r = RStr[PosR];
			char a = AStr[PosA];
			if (toupper(r) == toupper(a))
				++Ids;
			++PosR;
			++PosA;
			break;
			}

		case 'D':
			++PosR;
			break;

		case 'I':
			++PosA;
			break;

		default:
			asserta(false);
			}
		}
	asserta(PosR == LR);
	asserta(PosA == LA);
	double PctId = GetPct(Ids, PairColCount);

	WriteAlnPretty(g_fLog, (const byte *) RStr.c_str(), (const byte *) AStr.c_str(),
	  PathStr.c_str());
	ProgressLog("ref %s, add %s (%.1f%% id)\n", RefLabel.c_str(), AddLabel.c_str(), PctId);

	string PathStrXYB;
	for (uint i = 0; i < PairColCount; ++i)
		{
		char c = PathStr[i];
		if (c == 'M')
			PathStrXYB.push_back('B');
		else if (c == 'D')
			PathStrXYB.push_back('Y');
		else if (c == 'I')
			PathStrXYB.push_back('X');
		else
			asserta(false);
		}

	vector<uint> FreshIndexToMSAIndex;
	FreshIndexToMSAIndex.push_back(RefSeqIndex);
	vector<string> PWPaths;
	PWPaths.push_back(PathStrXYB);

	TransAln TA;
	TA.Init(RefAln, SeqsToAdd, FreshIndexToMSAIndex, PWPaths);
	TA.MakeExtendedMSA();
	TA.m_ExtendedMSA->ToFasta(opt(output));
	ProgressLog("Done.\n");
	}
