#include "muscle.h"

void cmd_relabel()
	{
	MultiSequence M;
	M.FromFASTA(opt(relabel));

	FILE *f = OpenStdioFile(opt(labels2));
	string Line;
	vector<string> Fields;
	map<string, string> OldLabelToNewLabel;
	uint LabelCount = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 2)
			Die("Expected 2 fields in line '%s'", Line.c_str());
		const string &OldLabel = Fields[0];
		const string &NewLabel = Fields[1];
		if (OldLabelToNewLabel.find(OldLabel) != OldLabelToNewLabel.end())
			Die("Dupe label >%s", OldLabel.c_str());
		OldLabelToNewLabel[OldLabel] = NewLabel;
		++LabelCount;
		}
	CloseStdioFile(f);

	const uint SeqCount = M.GetSeqCount();
	uint NotFound = 0;
	uint Found = 0;
	FILE *fOut = CreateStdioFile(opt(output));
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const Sequence *S = M.GetSequence(SeqIndex);
		const string &OldLabel = M.GetLabel(SeqIndex);
		map<string, string>::const_iterator p =
		  OldLabelToNewLabel.find(OldLabel);
		if (p == OldLabelToNewLabel.end())
			{
			S->WriteMFA(fOut);
			++NotFound;
			if (NotFound < 10)
				ProgressLog("Not found >%s\n", OldLabel.c_str());
			else if (NotFound == 10)
				ProgressLog("10+ Not found\n");
			continue;
			}
		const string &NewLabel = p->second;
		const byte *ByteSeq = S->GetBytePtr();
		uint L = S->GetLength();
		SeqToFasta(fOut, ByteSeq, L, NewLabel.c_str());
		}
	}
