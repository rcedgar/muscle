#include "muscle.h"
#include "ensemble.h"

void cmd_efa_bestconf()
	{
	const string &FileName = opt(efa_bestconf);

	Ensemble E;
	E.FromFile(FileName);

	const uint SeqCount = E.GetSeqCount();
	const uint MSACount = E.GetMSACount();
	const uint IxCount = E.GetIxCount();
	double AvgCols = double(IxCount)/MSACount;
	ProgressLog("%u seqs, %u MSAs, avg cols %.1f\n",
	  SeqCount, MSACount, AvgCols);

	ProgressLog("  MSA     Cols     N1   N1f  TotConf  MedConf  Name\n");
	//           12345  1234567  12345  1234  1234567  1234567
	uint BestMSAIndex_Total = 0;
	uint BestMSAIndex_Median = 0;
	string BestMSAName_Total;
	string BestMSAName_Median;
	double BestConf_Total = -1;
	double BestConf_Median = -1;
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *E.m_MSAs[MSAIndex];
		const string &Name = E.m_MSANames[MSAIndex];
		uint N1 = E.GetN1(MSAIndex);
		uint ColCount = M.GetColCount();
		double TotalConf = E.GetTotalConf(MSAIndex);
		double MedianConf = E.GetMedianConf(MSAIndex);
		if (TotalConf > BestConf_Total)
			{
			BestConf_Total = TotalConf;
			BestMSAIndex_Total = MSAIndex;
			BestMSAName_Total = Name;
			}
		if (MedianConf > BestConf_Median)
			{
			BestConf_Median = MedianConf;
			BestMSAIndex_Median = MSAIndex;
			BestMSAName_Median = Name;
			}
		double N1f = double(N1)/ColCount;
		ProgressLog("%5u  %7u  %5u  %4.2f  %7.3f  %7.4f  %s\n",
		  MSAIndex+1, ColCount, N1, N1f, TotalConf, MedianConf, Name.c_str());
		}

	ProgressLog("Best MSA, total  %u (%s)\n",
	  BestMSAIndex_Total+1, BestMSAName_Total.c_str());
	ProgressLog("Best MSA, median %u (%s)\n",
	  BestMSAIndex_Median+1, BestMSAName_Median.c_str());

	E.m_MSAs[BestMSAIndex_Median]->ToFASTAFile(opt(output));
	}
