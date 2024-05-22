#include "muscle.h"
#include "ensemble.h"
#include "qscorer.h"

void cmd_qscore_efa()
	{
	const string EfaFileName = opt(qscore_efa);
	const string RefFileName = opt(ref);
	const string OutputFileName = opt(output);
	double MaxGapFract = optd(max_gap_fract, 1.0);

	Ensemble E;
	E.FromFile(EfaFileName);

	MSA RefMSA;
	RefMSA.FromFASTAFile_PreserveCase(RefFileName);

	string RefName;
	GetBaseName(RefFileName, RefName);

	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;

	const uint MSACount = E.GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &TestMSA = E.GetMSA(MSAIndex);
		const string &TestName = E.GetMSAName(MSAIndex);
		QS.Run(TestName, TestMSA, RefMSA);
		ProgressLog("%s %s Q=%.4f TC=%.4f\n",
		  RefName.c_str(), TestName.c_str(), QS.m_Q, QS.m_TC);
		}
	}
