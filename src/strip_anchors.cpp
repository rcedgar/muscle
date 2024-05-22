#include "muscle.h"

void cmd_strip_anchors()
	{
	const string &InputFileName = opt(strip_anchors);
	MSA Aln;
	Progress("Reading %s ...", InputFileName.c_str());
	Aln.FromFASTAFile(InputFileName);
	Progress("done.\n");

	FILE *fOut = CreateStdioFile(opt(output));

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	asserta(SeqCount > 0);
	asserta(ColCount > 0);

	const char *Seq0 = Aln.GetSeqCharPtr(0);

	string AnchorStr;
	uint AnchorCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Seq0[Col];
		if (isdigit(c))
			{
			AnchorStr.push_back(c);
			++AnchorCount;
			}
		else
			AnchorStr.push_back(' ');
		}
	ProgressLog("%u anchors\n", AnchorCount);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Converting");
		Pf(fOut, ">%s\n", Aln.GetSeqName(SeqIndex));

		string SeqStr;
		const char *SeqCharPtr = Aln.GetSeqCharPtr(SeqIndex);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char a = AnchorStr[Col];
			char c = SeqCharPtr[Col];
			if (isdigit(a))
				{
				if (c != a)
					{
					string Label;
					Aln.GetSeqLabel(SeqIndex, Label);
					Log("\n");
					Log(">%s\n", Label.c_str());
					Log("AnchorStr %s\n", AnchorStr.c_str());
					Log("      Seq %s\n", SeqCharPtr);
					Die("Mis-aligned anchor");
					}
				}
			else
				SeqStr += c;
			}

		Pf(fOut, "%s\n", SeqStr.c_str());
		}

	CloseStdioFile(fOut);
	}
