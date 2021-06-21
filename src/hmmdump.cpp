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
		float Prob = EXP(Score);
		InitProbs.push_back(Prob);
		InitScores.push_back(Score);
		}

	const float InitProb_M = InitProbs[HMMSTATE_M];
	const float InitProb_ISX = InitProbs[HMMSTATE_ISX];
	const float InitProb_ISY = InitProbs[HMMSTATE_ISY];
	const float InitProb_ILX = InitProbs[HMMSTATE_ILX];
	const float InitProb_ILY = InitProbs[HMMSTATE_ILY];
	const float InitSum = InitProb_M +
	  InitProb_ISX + InitProb_ISY +
	  InitProb_ILX + InitProb_ILY;

	const float InitProb_IS = InitProb_ISX;
	const float InitProb_IL = InitProb_ILX;

	asserta(feq(InitProb_ISX, InitProb_ISY));
	asserta(feq(InitProb_ILX, InitProb_ILY));
	asserta(feq(InitSum, 1.0));

	fprintf(fOut, "\n");
	fprintf(fOut, "// Probs\n");
	fprintf(fOut, "const float InitProb_IM = %.5g;\n", InitProb_M);
	fprintf(fOut, "const float InitProb_IS = %.5g;\n", InitProb_ISX);
	fprintf(fOut, "const float InitProb_IL = %.5g;\n", InitProb_ILX);

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
			float Prob = EXP(Score);
			Sum += Prob;
			TransProbs[i][j] = Prob;
			TransScores[i][j] = Score;
			}
		asserta(feq(Sum, 1.0));
		}

// No transitions between different insert states
	asserta(TransProbs[HMMSTATE_ISX][HMMSTATE_ISY] == 0);
	asserta(TransProbs[HMMSTATE_ISX][HMMSTATE_ILY] == 0);
	asserta(TransProbs[HMMSTATE_ILX][HMMSTATE_ISY] == 0);
	asserta(TransProbs[HMMSTATE_ILX][HMMSTATE_ILY] == 0);

	asserta(TransProbs[HMMSTATE_ISY][HMMSTATE_ISX] == 0);
	asserta(TransProbs[HMMSTATE_ISY][HMMSTATE_ILX] == 0);
	asserta(TransProbs[HMMSTATE_ILY][HMMSTATE_ISX] == 0);
	asserta(TransProbs[HMMSTATE_ILY][HMMSTATE_ILX] == 0);

	asserta(TransProbs[HMMSTATE_M][HMMSTATE_ISX] ==
	  TransProbs[HMMSTATE_M][HMMSTATE_ISY]);

	asserta(TransProbs[HMMSTATE_M][HMMSTATE_ILX] ==
	  TransProbs[HMMSTATE_M][HMMSTATE_ILY]);

	asserta(TransProbs[HMMSTATE_ISX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_ISY][HMMSTATE_M]);

	asserta(TransProbs[HMMSTATE_ILX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_ILY][HMMSTATE_M]);

	asserta(TransProbs[HMMSTATE_ISX][HMMSTATE_M] ==
	  TransProbs[HMMSTATE_ISY][HMMSTATE_M]);

	const float TransProb_M_M = TransProbs[HMMSTATE_M][HMMSTATE_M];
	const float TransScore_M_M = TransScores[HMMSTATE_M][HMMSTATE_M];

	const float TransProb_M_IS = TransProbs[HMMSTATE_M][HMMSTATE_ISX];
	const float TransScore_M_IS = TransScores[HMMSTATE_M][HMMSTATE_ISX];

	const float TransProb_M_IL = TransProbs[HMMSTATE_M][HMMSTATE_ILX];
	const float TransScore_M_IL = TransScores[HMMSTATE_M][HMMSTATE_ILX];

	const float TransProb_IS_IS = TransProbs[HMMSTATE_ISX][HMMSTATE_ISX];
	const float TransScore_IS_IS = TransScores[HMMSTATE_ISX][HMMSTATE_ISX];

	const float TransProb_IL_IL = TransProbs[HMMSTATE_ILX][HMMSTATE_ILX];
	const float TransScore_IL_IL = TransScores[HMMSTATE_ILX][HMMSTATE_ILX];

	const float TransProb_IS_M = TransProbs[HMMSTATE_ISX][HMMSTATE_M];
	const float TransScore_IS_M = TransScores[HMMSTATE_ISX][HMMSTATE_M];

	const float TransProb_IL_M = TransProbs[HMMSTATE_ILX][HMMSTATE_M];
	const float TransScore_IL_M = TransScores[HMMSTATE_ILX][HMMSTATE_M];

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
		float Prob = EXP(Score);
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
			float Prob = EXP(Score);
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
	fprintf(fOut, "const float InitScore_IS = %.5g;\n", InitScores[HMMSTATE_ISX]);
	fprintf(fOut, "const float InitScore_IL = %.5g;\n", InitScores[HMMSTATE_ILX]);

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

	InitProbcons();
	
	PairHMM::WriteParamsReport(OutDir + "params_report.txt");

	HMMParams HP;
	HP.FromDefaults();
	HP.ToFile(OutDir + "hmm.tsv");

	HP.ToPairHMM();
	PairHMM::WriteParamsReport(OutDir + "params_report2.txt");

	HP.ToFile("hmm2.tsv");
	HP.FromFile("hmm2.tsv");
	HP.ToFile("hmm3.tsv");
	}
