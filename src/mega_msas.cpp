#include "muscle.h"

static char GetFeatureChar(bool IsAA, uint Letter)
	{
	if (IsAA)
		return g_CharToLetterAmino[Letter];
	if (Letter < 26)
		return 'A' + Letter;
	else if (Letter < 26*2)
		return 'a' + Letter;
	asserta(false);
	return '?';
	}

void cmd_mega_msas()
	{
	asserta(optset_input);
	const string &FaFN = opt(input);
	MSA Aln;
	Aln.FromFASTAFile_PreserveCase(FaFN);
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();

	MultiSequence InputSeqs;
	LoadInput(InputSeqs);
	const uint FeatureCount = Mega::GetFeatureCount();
	if (FeatureCount == 0)
		Die("No features in %s", g_Arg1.c_str());

	const string &OutputPrefix = opt(output);

	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		const string &FeatureName = Mega::m_FeatureNames[FeatureIdx];
		bool IsAA = (FeatureName == "AA");
		string OutputFN = OutputPrefix + FeatureName;

		FILE *f = CreateStdioFile(OutputFN);
		
		for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
			{
			const char *AARow = Aln.GetSeqCharPtr(SeqIdx);
			const string &Label = Aln.GetLabel(SeqIdx);
			const vector<vector<byte> > &Profile =
			  *Mega::GetProfileByLabel(Label);
			uint Pos = 0;
			string FeatureRow;
			for (uint Col = 0; Col < ColCount; ++Col)
				{
				char c = AARow[Col];
				if (isgap(c))
					FeatureRow += c;
				else
					{
					asserta(SIZE(Profile[Pos]) == FeatureCount);
					byte FeatureLetter = Profile[Pos][FeatureIdx];
					FeatureRow += GetFeatureChar(IsAA, FeatureLetter);
					++Pos;
					}
				}
			SeqToFasta(f, FeatureRow, Label);
			}

		CloseStdioFile(f);
		}
	}
