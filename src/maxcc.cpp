#include "muscle.h"
#include "ensemble.h"

void cmd_maxcc()
	{
	const string &InputFileName = opt(maxcc);
	const string &OutputFileName = opt(output);
	if (OutputFileName.empty())
		Die("Must set -output");

	Ensemble E;
	E.FromFile(InputFileName);
	const uint MSACount = E.GetMSACount();
	if (MSACount == 0)
		Die("Ensemble is empty");

	uint BestMSAIndex = 0;
	double BestConf = 0;
	double SumConf = 0;
	double MinConf = 0;
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &Name = E.GetMSAName(MSAIndex);
		const MSA &M = *E.m_MSAs[MSAIndex];
		double TotalConf = E.GetTotalConf(MSAIndex);
		ProgressStep(MSAIndex, MSACount, "%s (%.3g)", Name.c_str(), TotalConf);
		Log("%u/%u %s (%.3g)\n", MSAIndex, MSACount, Name.c_str(), TotalConf);
		if (MSAIndex == 0 || TotalConf < MinConf)
			MinConf = TotalConf;
		SumConf += TotalConf;
		if (TotalConf >= BestConf)
			{
			BestMSAIndex = MSAIndex;
			BestConf = TotalConf;
			}
		}

	const string &BestName = E.GetMSAName(BestMSAIndex);
	const MSA *BestMSA = E.m_MSAs[BestMSAIndex];
	double AvgConf = SumConf/MSACount;
	ProgressLog("CC min %.3g, avg %.3g, max %.3g, best %s\n",
	  MinConf, AvgConf, BestConf, BestName.c_str());

	asserta(BestMSA != 0);
	BestMSA->ToFASTAFile(OutputFileName);
	}
