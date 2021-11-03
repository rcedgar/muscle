#include "muscle.h"

static char GetConsChar(const MultiSequence &MSA, uint ColIndex)
	{
	asserta(g_AlphaSize == 4 || g_AlphaSize == 20);
	vector<uint> Counts(g_AlphaSize+1);
	const uint ColCount = MSA.GetColCount();
	const uint SeqCount = MSA.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = MSA.GetChar(SeqIndex, ColIndex);
		if (isgap(c))
			{
			++(Counts[g_AlphaSize]);
			continue;
			}
		uint Letter = CharToLetter(c);
		if (Letter < g_AlphaSize)
			++(Counts[Letter]);
		}

	uint MaxCount = 0;
	uint MaxLetter = 0;
	for (uint Letter = 0; Letter <= g_AlphaSize; ++Letter)
		{
		uint Count = Counts[Letter];
		if (Count > MaxCount)
			{
			MaxCount = Count;
			MaxLetter = Letter;
			}
		}
	if (MaxLetter == g_AlphaSize)
		return '-';
	char ConsChar = LetterToChar(MaxLetter);
	return ConsChar;
	}

void GetConsensusSequence(const MultiSequence &MSA, string &Seq)
	{
	Seq.clear();

	const uint SeqCount = MSA.GetSeqCount();
	const uint ColCount = MSA.GetColCount();

	vector<uint> Freqs(g_AlphaSize);
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		char c = GetConsChar(MSA, ColIndex);
		if (c != '-')
			Seq += c;
		}
	}

void cmd_consseq()
	{
	const string &MSAFileName = opt(consseq);
	const string &OutputFileName = opt(output);
	string Label = "CONSENSUS";
	if (optset_label)
		Label = opt(label);

	MultiSequence MSA;
	MSA.FromFASTA(MSAFileName);

	string ConsSeq;
	GetConsensusSequence(MSA, ConsSeq);

	FILE *fOut = CreateStdioFile(OutputFileName);
	SeqToFasta(fOut, ConsSeq, Label);
	CloseStdioFile(fOut);
	}
