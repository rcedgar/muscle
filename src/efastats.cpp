#include "muscle.h"
#include "ensemble.h"
#include "qscorer.h"

static void CmpRef(const Ensemble &E, const MSA &RefMSA,
  double MaxGapFract, vector<double> &Qs, vector<double> &TCs)
	{
	Qs.clear();
	TCs.clear();

	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;
	const uint MSACount = E.GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &TestMSA = E.GetMSA(MSAIndex);
		const string &TestName = E.GetMSAName(MSAIndex);
		QS.Run(TestName, TestMSA, RefMSA);
		double Q = QS.m_Q;
		double TC = QS.m_TC;
		Qs.push_back(Q);
		TCs.push_back(TC);
		}
	}

void cmd_efastats()
	{
	const string &InputFileName = opt(efastats);
	double MaxGapFract = optd(max_gap_fract, 0.5);
	const string &RefFileName = opt(ref);

	Ensemble E;
	E.FromFile(InputFileName);

	vector<double> Qs;
	vector<double> TCs;
	if (optset_ref)
		{
		MSA RefMSA;
		RefMSA.FromFASTAFile(opt(ref));
		CmpRef(E, RefMSA, MaxGapFract, Qs, TCs);
		}

	const uint SeqCount = E.GetSeqCount();
	const uint MSACount = E.GetMSACount();
	const uint IxCount = E.GetIxCount();

	double D_LetterPairs;
	double D_Columns;
	E.GetDispersion(MaxGapFract, D_LetterPairs, D_Columns);

	vector<double> CCs;
	double AvgCols = double(IxCount)/MSACount;

	ProgressLog("  MSA     Cols     N1   N1f  Conf     CC");
	//           12345  1234567  12345  1234  1234  12345
	if (optset_ref)
		ProgressLog("       Q      TC");
		//             123456  123456
	ProgressLog("  Name\n");
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *E.m_MSAs[MSAIndex];
		const string &Name = E.m_MSANames[MSAIndex];
		uint N1 = E.GetN1(MSAIndex);
		uint ColCount = M.GetColCount();
		double TotalConf = E.GetTotalConf(MSAIndex);
		double N1f = double(N1)/ColCount;
		double CC = TotalConf/ColCount;
		CCs.push_back(CC);
		ProgressLog("%5u  %7u  %5u  %4.2f  %4.2f  %5.3f",
		  MSAIndex+1, ColCount, N1, N1f, TotalConf, CC);
		if (optset_ref)
			ProgressLog("  %6.4f  %6.4f", Qs[MSAIndex], TCs[MSAIndex]);
		ProgressLog("  %s\n", Name.c_str());
		}

	sort(CCs.begin(), CCs.end());
	double MedianCC = CCs[MSACount/2];

	Progress("%u seqs, %u MSAs, avg cols %.1f, D_LP %.3g, D_Cols %.3g, CC %.3g",
	  SeqCount, MSACount, AvgCols, D_LetterPairs, D_Columns, MedianCC);
	Log("@SUMMARY input=%s D_LP=%.4f D_Cols=%.4f CC=%.4f",
	  InputFileName.c_str(), D_LetterPairs, D_Columns, MedianCC);
	if (optset_ref)
		{
		sort(Qs.begin(), Qs.end());
		sort(TCs.begin(), TCs.end());
		asserta(SIZE(Qs) == MSACount);
		asserta(SIZE(TCs) == MSACount);
		double MedianQ = Qs[MSACount/2];
		double MedianTC = TCs[MSACount/2];
		double E_LP = 1 - MedianQ;
		double E_Cols = 1 - MedianTC;
		Progress(" E_LP %.4f, E_Cols %.4f", E_LP, E_Cols);
		Log(" E_LP=%.4f E_Cols=%.4f", E_LP, E_Cols);
		}
	Progress("\n");
	Log("\n");
	}
