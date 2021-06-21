#include "myutils.h"
#include "probcons.h"

void cmd_pp()
	{
	const string &MSAFileName1 = opt(pp);
	const string &MSAFileName2 = opt(input2);

	int TargetPairCount = 1024;
	if (optset_paircount)
		TargetPairCount = int(opt(paircount));

	InitProbcons();

	string BaseName1;
	string BaseName2;
	GetBaseName(MSAFileName1, BaseName1);
	GetBaseName(MSAFileName2, BaseName2);

	MultiSequence MSA1;
	MultiSequence MSA2;
	Progress("Reading %s ...", BaseName1.c_str());
	MSA1.LoadMFA(MSAFileName1, false);
	Progress("done.\n");

	Progress("Reading %s ...", BaseName2.c_str());
	MSA2.LoadMFA(MSAFileName2, false);
	Progress("done.\n");

	vector<char> Path;
	string ProgressStr = BaseName1 + "+" + BaseName2;
	float EA = 
		AlignMSAs(ProgressStr, MSA1, MSA2, TargetPairCount, Path);
	ProgressLog("%s: EA %.4g\n", ProgressStr.c_str(), EA);

	if (optset_output)
		{
		MultiSequence MSA12;

		for (int SeqIndex = 0; SeqIndex < MSA1.GetNumSequences(); ++SeqIndex)
			{
			const Sequence *Seq1 = MSA1.GetSequence(SeqIndex);
			Sequence *AlignedSeq1 = Seq1->AddGaps(&Path, 'X');
			MSA12.AddSequence(AlignedSeq1);
			}

		for (int SeqIndex = 0; SeqIndex < MSA2.GetNumSequences(); ++SeqIndex)
			{
			const Sequence *Seq2 = MSA2.GetSequence(SeqIndex);
			Sequence *AlignedSeq2 = Seq2->AddGaps(&Path, 'Y');
			MSA12.AddSequence(AlignedSeq2);
			}

		MSA12.WriteMFA(opt(output));
		}
	}
