#include "muscle.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path);

void cmd_transalnref()
	{
	asserta(optset_input2);
	asserta(optset_output);
	asserta(optset_label);

	const string &RefAlnFN = g_Arg1;
	const string &AddFastaFN = opt(input2);
	const string &RefLabel = opt(label);

	MSA RefAln;
	RefAln.FromFASTAFile(RefAlnFN);
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

	string R;
	RefAln.GetUngappedSeqStr(RefSeqIndex, R);
	const uint LR = SIZE(R);

	MSA SeqsToAdd;
	SeqsToAdd.FromFASTAFile(AddFastaFN);
	if (SeqsToAdd.GetSeqCount() != 1)
		Die("-input2 must have exactly one sequence");

	string A;
	SeqsToAdd.GetUngappedSeqStr(0, A);
	const uint LA = SIZE(A);

	PathInfo *PI = ObjMgr::GetPathInfo();
	XDPMem Mem;
	ViterbiFastMem(Mem, (const byte *) R.c_str(), LR, (const byte *) A.c_str(), LA, *PI);

	string PathStr;
	PI->GetPathStr(PathStr);
	Log("Path=%s\n", PathStr.c_str());
	WriteAlnPretty(g_fLog, (const byte *) R.c_str(), (const byte *) A.c_str(),
	  PathStr.c_str());
	}
