#include "muscle.h"
#include "hmmparams.h"

void PairHMM::WriteParamsReport(const string &FileName)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	WriteParamsReport(f);
	CloseStdioFile(f);
	}

void PairHMM::WriteParamsReport(FILE *fOut)
	{
	if (fOut == 0)
		return;

	vector<float> InitProbs;
	vector<float> InitScores;
	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		{
		float Score = PairHMM::m_StartScore[i];
		float Prob = exp(Score);
		InitProbs.push_back(Prob);
		InitScores.push_back(Score);
		}

	const float InitProb_M = InitProbs[HMMSTATE_M];
	const float InitProb_IX = InitProbs[HMMSTATE_IX];
	const float InitProb_IY = InitProbs[HMMSTATE_IY];
	const float InitProb_JX = InitProbs[HMMSTATE_JX];
	const float InitProb_JY = InitProbs[HMMSTATE_JY];
	const float InitSum = InitProb_M +
	  InitProb_IX + InitProb_IY +
	  InitProb_JX + InitProb_JY;

	const float InitProb_IS = InitProb_IX;
	const float InitProb_IL = InitProb_JX;

	asserta(feq(InitProb_IX, InitProb_IY));
	asserta(feq(InitProb_JX, InitProb_JY));
	asserta(feq(InitSum, 1.0));

	fprintf(fOut, "\n");
	fprintf(fOut, "// Probs\n");
	fprintf(fOut, "const float InitProb_IM = %.5g;\n", InitProb_M);
	fprintf(fOut, "const float InitProb_IS = %.5g;\n", InitProb_IX);
	fprintf(fOut, "const float InitProb_IL = %.5g;\n", InitProb_JX);

	vector<vector<float> > TransProbs(HMMSTATE_COUNT);
	vector<vector<float> > TransScores(HMMSTATE_COUNT);
	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		{
		TransProbs[i].resize(HMMSTATE_COUNT);
		TransScores[i].resize(HMMSTATE_COUNT);

		float Sum = 0;
		for (uint j = 0; j < HMMSTATE_COUNT; ++j)
			{
			float Score = PairHMM::m_TransScore[i][j];
			float Prob = exp(Score);
			Sum += Prob;
			TransProbs[i][j] = Prob;
			TransScores[i][j] = Score;
			}
		asserta(feq(Sum, 1.0));
		}

// No transitions between different insert states
	asserta(TransProbs[HMMSTATE_IX][HMMSTATE_IY] == 0);
	asserta(TransProbs[HMMSTATE_IX][HMMSTATE_JY] == 0);
	asserta(TransProbs[HMMSTATE_JX][HMMSTATE_IY] == 0);
	asserta(TransProbs[HMMSTATE_JX][HMMSTATE_JY] == 0);

	asserta(TransProbs[HMMSTATE_IY][HMMSTATE_IX] == 0);
	asserta(TransProbs[HMMSTATE_IY][HMMSTATE_JX] == 0);
	asserta(TransProbs[HMMSTATE_JY][HMMSTATE_IX] == 0);
	asserta(TransProbs[HMMSTATE_JY][HMMSTATE_JX] == 0);

	asserta(TransProbs[HMMSTATE_M][HMMSTATE_IX] ==
	  TransProbs[HMMSTATE_M][HMMSTATE_IY]);

	asserta(TransProbs[HMMSTATE_M][HMMSTATE_JX] ==
	  TransProbs[HMMSTATE_M][HMMSTATE_JY]);

	asserta(TransProbs[HMMSTATE_IX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_IY][HMMSTATE_M]);

	asserta(TransProbs[HMMSTATE_JX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_JY][HMMSTATE_M]);

	asserta(TransProbs[HMMSTATE_IX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_IY][HMMSTATE_M]);

	const float TransProb_M_M = TransProbs[HMMSTATE_M][HMMSTATE_M];
	const float TransScore_M_M = TransScores[HMMSTATE_M][HMMSTATE_M];

	const float TransProb_M_IS = TransProbs[HMMSTATE_M][HMMSTATE_IX];
	const float TransScore_M_IS = TransScores[HMMSTATE_M][HMMSTATE_IX];

	const float TransProb_M_IL = TransProbs[HMMSTATE_M][HMMSTATE_JX];
	const float TransScore_M_IL = TransScores[HMMSTATE_M][HMMSTATE_JX];

	const float TransProb_IS_IS = TransProbs[HMMSTATE_IX][HMMSTATE_IX];
	const float TransScore_IS_IS = TransScores[HMMSTATE_IX][HMMSTATE_IX];

	const float TransProb_IL_IL = TransProbs[HMMSTATE_JX][HMMSTATE_JX];
	const float TransScore_IL_IL = TransScores[HMMSTATE_JX][HMMSTATE_JX];

	const float TransProb_IS_M = TransProbs[HMMSTATE_IX][HMMSTATE_M];
	const float TransScore_IS_M = TransScores[HMMSTATE_IX][HMMSTATE_M];

	const float TransProb_IL_M = TransProbs[HMMSTATE_JX][HMMSTATE_M];
	const float TransScore_IL_M = TransScores[HMMSTATE_JX][HMMSTATE_M];

	asserta(feq(InitProb_M + 2*InitProb_IS + 2*InitProb_IL, 1.0));
	asserta(feq(TransProb_IS_IS + TransProb_IS_M, 1.0));
	asserta(feq(TransProb_IL_IL + TransProb_IL_M, 1.0));
	asserta(feq(TransProb_M_M + 2*TransProb_M_IS + 2*TransProb_M_IL, 1.0));

	fprintf(fOut, "\n");
	fprintf(fOut, "const float TransProb_M_M   = %.5g;\n", TransProb_M_M);
	fprintf(fOut, "const float TransProb_M_IS  = %.5g;\n", TransProb_M_IS);
	fprintf(fOut, "const float TransProb_M_IL  = %.5g;\n", TransProb_M_IL);
	fprintf(fOut, "const float TransProb_IS_IS = %.5g;\n", TransProb_IS_IS);
	fprintf(fOut, "const float TransProb_IL_IL = %.5g;\n", TransProb_IL_IL);
	fprintf(fOut, "const float TransProb_IS_M  = %.5g;\n", TransProb_IS_M);
	fprintf(fOut, "const float TransProb_IL_M  = %.5g;\n", TransProb_IL_M);

	const string A = "ACDEFGHIKLMNPQRSTVWY";
	asserta(SIZE(A) == 20);
	vector<float> InsProbs(20);
	vector<float> InsScores(20);
	float Sum = 0;
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		uint Letter = uint(a);
		float Score = PairHMM::m_InsScore[Letter];
		float Prob = exp(Score);
		InsProbs[i] = Prob;
		InsScores[i] = Score;
		Sum += Prob;
		}
	asserta(feq(Sum, 1.0));

	vector<vector<float> > EmitProbs(20);
	vector<vector<float> > EmitScores(20);
	Sum = 0;
	for (uint i = 0; i < 20; ++i)
		{
		EmitProbs[i].resize(20);
		EmitScores[i].resize(20);

		char a = A[i];
		uint Letter_a = uint(a);
		for (uint j = 0; j < 20; ++j)
			{
			char b = A[j];
			uint Letter_b = uint(b);
			float Score = PairHMM::m_MatchScore[Letter_a][Letter_b];
			float Prob = exp(Score);
			EmitProbs[i][j] = Prob;
			EmitScores[i][j] = Score;
			Sum += Prob;
			}
		}
	asserta(feq(Sum, 1.0));

	fprintf(fOut, "\n");
	fprintf(fOut, "const float InsProbs[20] =\n");
	fprintf(fOut, "	{\n");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		float Prob = InsProbs[i];
		fprintf(fOut, "	%.5g,	// %c\n", Prob, a);
		}
	fprintf(fOut, "	};\n");

	fprintf(fOut, "\n");
	fprintf(fOut, "const float EmitProbs[20][20] =\n");
	fprintf(fOut, "	{\n");
	fprintf(fOut, "//	      ");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		fprintf(fOut, "        %c", a);
		}
	fprintf(fOut, "\n");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		fprintf(fOut, "/* %c */ { ", a);
		for (uint j = 0; j < 20; ++j)
			{
			float Prob = EmitProbs[i][j];
			fprintf(fOut, " %.5g", Prob);
			}
		fprintf(fOut, " } // %c\n", a);
		}
	fprintf(fOut, "	};\n");

///////////////////////////////////////////////////////////////////
// Scores
///////////////////////////////////////////////////////////////////
	fprintf(fOut, "\n");
	fprintf(fOut, "// Scores\n");
	fprintf(fOut, "const float InitScore_IM = %.5g;\n", InitScores[HMMSTATE_M]);
	fprintf(fOut, "const float InitScore_IS = %.5g;\n", InitScores[HMMSTATE_IX]);
	fprintf(fOut, "const float InitScore_IL = %.5g;\n", InitScores[HMMSTATE_JX]);

	fprintf(fOut, "\n");
	fprintf(fOut, "const float TransScore_M_M   = %.5g;\n", TransScore_M_M);
	fprintf(fOut, "const float TransScore_M_IS  = %.5g;\n", TransScore_M_IS);
	fprintf(fOut, "const float TransScore_M_IL  = %.5g;\n", TransScore_M_IL);
	fprintf(fOut, "const float TransScore_IS_IS = %.5g;\n", TransScore_IS_IS);
	fprintf(fOut, "const float TransScore_IL_IL = %.5g;\n", TransScore_IL_IL);
	fprintf(fOut, "const float TransScore_IS_M  = %.5g;\n", TransScore_IS_M);
	fprintf(fOut, "const float TransScore_IL_M  = %.5g;\n", TransScore_IL_M);

	fprintf(fOut, "\n");
	fprintf(fOut, "const float InsScores[20] =\n");
	fprintf(fOut, "	{\n");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		float Score = InsScores[i];
		fprintf(fOut, "	%.5g,	// %c\n", Score, a);
		}
	fprintf(fOut, "	};\n");

	fprintf(fOut, "\n");
	fprintf(fOut, "const float EmitScores[20][20] =\n");
	fprintf(fOut, "	{\n");
	fprintf(fOut, "//	      ");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		fprintf(fOut, "        %c", a);
		}
	fprintf(fOut, "\n");
	for (uint i = 0; i < 20; ++i)
		{
		char a = A[i];
		fprintf(fOut, "/* %c */ { ", a);
		for (uint j = 0; j < 20; ++j)
			{
			float Prob = EmitScores[i][j];
			fprintf(fOut, " %.5g", Prob);
			}
		fprintf(fOut, " } // %c\n", a);
		}
	fprintf(fOut, "	};\n");
	}

void cmd_hmmdump()
	{
	string OutDir = opt(hmmdump);
	Dirize(OutDir);

	SetAlpha(ALPHA_Amino);
	InitProbcons();
	
	PairHMM::WriteParamsReport(OutDir + "params_report.txt");

	bool Nucleo = opt(nt);

	HMMParams HP;
	HP.FromDefaults(Nucleo);
	HP.ToFile(OutDir + "hmm.tsv");

	HP.CmdLineUpdate();
	HP.ToPairHMM();
	PairHMM::WriteParamsReport(OutDir + "params_report2.txt");

	HP.ToFile(OutDir + "hmm2.tsv");
	HP.FromFile(OutDir + "hmm2.tsv");
	HP.ToFile(OutDir + "hmm3.tsv");

	HMMParams SA;
	HP.ToSingleAffineProbs(SA);
	SA.ToFile(OutDir + "sa.hmm");
	}
