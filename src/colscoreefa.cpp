#include "muscle.h"
#include "ensemble.h"

#pragma warning(3: 4365)	// signed/unsigned conversion

static const uint MAXBIN = 10;

static uint GetBin(double Conf)
	{
	asserta(Conf > 0 && Conf <= 1);
	if (Conf == 1)
		return MAXBIN;
	uint Bin = uint(Conf*10);
	asserta(Bin >= 0 && Bin < MAXBIN);
	return Bin;
	}

void cmd_colscore_efa()
	{
	const string EfaFileName = opt(colscore_efa);
	const string RefFileName = opt(ref);
	const string OutputFileName = opt(output);
	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	Ensemble E;
	E.FromFile(EfaFileName);

	MSA Ref;
	Ref.FromFASTAFile_PreserveCase(RefFileName);

	FILE *fOut = CreateStdioFile(OutputFileName);

	const uint MSACount = E.GetMSACount();
	E.SortMSA(Ref);
	
	set<uint> RefUniqueIxs;
	E.GetRefUniqueIxs(Ref, RefUniqueIxs, MaxGapFract);
	const uint RefIxCount = SIZE(RefUniqueIxs);

	uint RefUpperColCount = 0;
	const uint RefColCount = Ref.GetColCount();
	for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
		if (Ref.ColIsUpper(RefColIndex, MaxGapFract))
			++RefUpperColCount;

	set<pair<uint, int> > RefPosSet;
	E.GetRefPosSet(Ref, MaxGapFract, RefPosSet);

	vector<uint> BinToCount(MAXBIN+1);
	vector<uint> BinToCorrectCount(MAXBIN+1);
	double SumTC = 0;
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		vector<uint> TestUniqueIxs;
		vector<double> Confs;
		E.GetTestUniqueIxs(MSAIndex, RefPosSet, TestUniqueIxs, Confs);

		const uint TestIxCount = SIZE(TestUniqueIxs);
		uint CorrectCount = 0;
		for (uint i = 0; i < TestIxCount; ++i)
			{
			uint TestUniqueIx = TestUniqueIxs[i];
			double Conf = Confs[i];
			uint Bin = GetBin(Conf);
			asserta(Bin <= MAXBIN);
			++BinToCount[Bin];
			bool Correct =
			  (RefUniqueIxs.find(TestUniqueIx) != RefUniqueIxs.end());
			if (Correct)
				{
				++CorrectCount;
				++BinToCorrectCount[Bin];
				}
//			Pf(fOut, "col	%c	%.4f\n", tof(Correct), Conf);
			}
		double TC = double(CorrectCount)/RefUpperColCount;
		SumTC += TC;
		//Pf(fOut, "tc	%u	%.4f\n", MSAIndex, TC);
		}
	double MeanTC = SumTC/MSACount;
	Pf(fOut, "meantc	%.4f\n", MeanTC);
	ProgressLog("Mean TC %.4f\n", MeanTC);

	ProgressLog("Bins ");
	for (uint Bin = 0; Bin <= MAXBIN; ++Bin)
		{
		uint Count = BinToCount[Bin];
		uint CorrectCount = BinToCorrectCount[Bin];
		asserta(CorrectCount <= Count);
		double P = 0;
		if (Count > 0)
			P = double(CorrectCount)/Count;
		Pf(fOut, "bin	%u	%u	%u	%.4f\n",
		  Bin, Count, CorrectCount, P);
		ProgressLog(" %.2f", P);
		}
	ProgressLog("\n");

	CloseStdioFile(fOut);
	}
