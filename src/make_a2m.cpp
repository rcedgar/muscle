#include "muscle.h"

extern bool g_FASTA_AllowDigits;

void cmd_make_a2m()
	{
	const string &InputFileName = opt(make_a2m);
	MSA Aln;
	Progress("Reading %s ...", InputFileName.c_str());
	Aln.FromFASTAFile(InputFileName);
	Progress("done.\n");

	FILE *fOut = CreateStdioFile(opt(output));

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	asserta(SeqCount > 0);
	asserta(ColCount > 0);

	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	vector<bool> MatchVec;
	uint MatchCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = Aln.GetGapCount(Col);
		double GapFract = double(GapCount)/double(SeqCount);
		bool IsMatch = (GapFract <= MaxGapFract);
		MatchVec.push_back(IsMatch);
		if (IsMatch)
			++MatchCount;
		}

	ProgressLog("%u match cols (%.1f%%), maxgapfract %.3f\n",
	  MatchCount, GetPct(MatchCount, ColCount), MaxGapFract);
	asserta(MatchCount > 0);

	asserta(SIZE(MatchVec) == ColCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Converting");
		Pf(fOut, ">%s\n", Aln.GetSeqName(SeqIndex));

		string SeqStr;
		const char *SeqCharPtr = Aln.GetSeqCharPtr(SeqIndex);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c = SeqCharPtr[Col];

			bool IsMatch = MatchVec[Col];
		
			if (IsMatch)
				{
				if (c == '.' || c == '-')
					c = '-';
				else if (isalpha(c))
					c = toupper(c);
				else if (g_FASTA_AllowDigits && isdigit(c))
					;
				else
					Die("Bad char 0x%02x", c);
				}
			else
				{
				asserta(!IsMatch);
				if (c == '.' || c == '-')
					continue;
				c = tolower(c);
				}
			SeqStr += c;
			}

		Pf(fOut, "%s\n", SeqStr.c_str());
		}

	CloseStdioFile(fOut);
	}
