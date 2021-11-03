#include "muscle.h"
#include "pprog.h"

void cmd_eadistmx_msas()
	{
	const string &FileName = opt(eadistmx_msas);
	vector<string> MSAFileNames;
	ReadStringsFromFile(FileName, MSAFileNames);
	const uint MSACount = SIZE(MSAFileNames);
	asserta(optset_output);
	FILE *f = CreateStdioFile(opt(output));

	PProg PP;
	if (optset_paircount)
		PP.m_TargetPairCount = opt(paircount);
	bool IsNucleo;
	PP.LoadMSAs(MSAFileNames, IsNucleo);
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);

	InitProbcons();
	PP.AlignAllInputPairs();

	vector<vector<float> > &ScoreMx = PP.m_ScoreMx;
	for (uint i = 0; i < MSACount; ++i)
		{
		asserta(i < SIZE(ScoreMx));
		const char *Labeli = PP.GetMSALabel(i).c_str();
		for (uint j = i+1; j < MSACount; ++j)
			{
			asserta(j < SIZE(ScoreMx[i]));
			const char *Labelj = PP.GetMSALabel(j).c_str();
			float Score = ScoreMx[i][j];
			fprintf(f, "%s\t%s\t%.4f\n", Labeli, Labelj, Score);
			}
		}

	CloseStdioFile(f);
	}
