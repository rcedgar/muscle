#include "muscle.h"
#include "qscorer.h"
#include "ensemble.h"

char ConfToChar_1(double Conf);
void WriteLetterConfHTML(const string &FileName,
  const MSA &Ref, const MSA &ConfAln);
void WriteLetterConfJalView(const string &FileName,
  const MSA &Ref, const MSA &ConfAln);

void Ensemble::GetLetterConfsVec(const MSA &Ref, double MaxGapFract,
  vector<vector<double> > &LetterConfsVec) const
	{
	QScorer QS;

	vector<vector<uint> > LetterCountsVec;
	const uint TestMSACount = GetMSACount();
	for (uint TestMSAIndex = 0; TestMSAIndex < TestMSACount; ++TestMSAIndex)
		{
		const MSA &Test = GetMSA(TestMSAIndex);
		QS.Run("Ensemble::GetLetterConfsVec()", Test, Ref);
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
			double LetterConf = 0;
			if (c == '-' || c == '.')
				asserta(n == 0);
			else
				LetterConf = double(n)/TestMSACount;
			LetterConfsVec[RefSeqIndex].push_back(LetterConf);
			}
		}
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

	vector<vector<double> > LetterConfsVec;
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
			char c = ConfAln.GetChar(RefSeqIndex, RefColIndex);
			if (isgap(c))
				continue;
			double Conf = LetterConfsVec[RefSeqIndex][RefColIndex];
			char ConfChar = ConfToChar_1(Conf);
			ConfAln.SetChar(RefSeqIndex, RefColIndex, ConfChar);
			}
		}
	ConfAln.ToFASTAFile(opt(output));
	WriteLetterConfHTML(opt(html), RefPC, ConfAln);
	WriteLetterConfJalView(opt(jalview), RefPC, ConfAln);
	}

void cmd_addletterconfseq()
	{
	const string EnsembleFileName = g_Arg1;
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

	vector<vector<double> > LetterConfsVec;
	E.GetLetterConfsVec(Ref, MaxGapFract, LetterConfsVec);
	asserta(SIZE(LetterConfsVec) == RefSeqCount);
	asserta(SIZE(LetterConfsVec[0]) == RefColCount);

	string ConfRow;
	for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
		{
		double SumConf = 0;
		for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
			{
			const string RefLabel = Ref.GetSeqName(RefSeqIndex);
			string SeqStr;
			SumConf += LetterConfsVec[RefSeqIndex][RefColIndex];
			}
		double MeanConf = SumConf/RefSeqCount;
		char ConfChar = ConfToChar_1(MeanConf);
		ConfRow += ConfChar;
		}

	FILE *fOut = CreateStdioFile(opt(output));
	SeqToFasta(fOut, ConfRow.c_str(), RefColCount, "_letterconf_");
	RefPC.ToFASTAFile(fOut);
	CloseStdioFile(fOut);
	}
