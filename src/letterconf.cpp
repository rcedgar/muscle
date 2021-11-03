#include "muscle.h"
#include "qscorer.h"
#include "ensemble.h"

void WriteLetterConfHTML(const string &FileName,
  const MSA &Ref, const MSA &ConfAln);
void WriteLetterConfJalView(const string &FileName,
  const MSA &Ref, const MSA &ConfAln);

void Ensemble::GetLetterConfsVec(const MSA &Ref, double MaxGapFract,
  vector<vector<uint> > &LetterConfsVec) const
	{
	QScorer QS;

	vector<vector<uint> > LetterCountsVec;
	const uint TestMSACount = GetMSACount();
	for (uint TestMSAIndex = 0; TestMSAIndex < TestMSACount; ++TestMSAIndex)
		{
		const MSA &Test = GetMSA(TestMSAIndex);
		QS.Run(Test, Ref);
		QS.UpdateRefLetterCounts(LetterCountsVec);
		}

	const uint RefSeqCount = QS.GetRefSeqCount();
	const uint RefColCount = QS.GetRefColCount();
	asserta(SIZE(LetterCountsVec) == RefSeqCount);
	asserta(SIZE(LetterCountsVec[0]) == RefColCount);

	LetterConfsVec.resize(RefSeqCount);
	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
			{
			uint n = LetterCountsVec[RefSeqIndex][RefColIndex];
			char c = Ref.GetChar(RefSeqIndex, RefColIndex);
			uint LetterConf = UINT_MAX;
			if (c == '-' || c == '.')
				asserta(n == 0);
			else
				LetterConf = (n*9)/TestMSACount;
			LetterConfsVec[RefSeqIndex].push_back(LetterConf);
			}
		}
	}

static char GetConfChar(uint n)
	{
	if (n == UINT_MAX)
		return '-';
	asserta(n <= 9);
	return '0' + n;
	}

void cmd_letterconf()
	{
	const string EnsembleFileName = opt(letterconf);
	const string RefFileName = opt(ref);

	double MaxGapFract = optd(max_gap_fract, 1.0);

	string Name;
	GetBaseName(RefFileName.c_str(), Name);

	MSA Ref;
	MSA RefPC;
	Ref.FromFASTAFile(RefFileName);
	RefPC.FromFASTAFile_PreserveCase(RefFileName);
	const uint RefSeqCount = Ref.GetSeqCount();
	const uint RefColCount = Ref.GetColCount();

	Ensemble E;
	E.FromFile(EnsembleFileName);

	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;

	vector<vector<uint> > LetterConfsVec;
	E.GetLetterConfsVec(Ref, MaxGapFract, LetterConfsVec);
	asserta(SIZE(LetterConfsVec) == RefSeqCount);
	asserta(SIZE(LetterConfsVec[0]) == RefColCount);

	MSA ConfAln;
	ConfAln.Copy(Ref);

	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const string RefLabel = Ref.GetSeqName(RefSeqIndex);
		string SeqStr;
		for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
			{
			uint n = LetterConfsVec[RefSeqIndex][RefColIndex];
			asserta(n <= 9 || n == UINT_MAX);
			char ConfChar = GetConfChar(n);
			ConfAln.SetChar(RefSeqIndex, RefColIndex, ConfChar);
			}
		}
	ConfAln.ToFASTAFile(opt(output));
	WriteLetterConfHTML(opt(html), RefPC, ConfAln);
	WriteLetterConfJalView(opt(jalview), RefPC, ConfAln);
	}
