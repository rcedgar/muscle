#include "muscle.h"
#include "ensemble.h"
#include "qscorer.h"

static char ConfToChar1(double Conf)
	{
	asserta(Conf >= 0 && Conf <= 1);
	uint Tenth = uint(Conf*10);
	asserta(Tenth >= 0 && Tenth <= 10);
	if (Tenth == 10)
		return '+';
	return '0' + Tenth;
	}

static char ConfToChar2(double Conf)
	{
	asserta(Conf >= 0 && Conf <= 1);
	uint H = uint(Conf*100);
	asserta(H >= 0 && H <= 100);
	if (H == 100)
		return '+';
	return '0' + H%10;
	}

char ConfToChar_1(double Conf)
	{
	asserta(Conf >= 0 && Conf <= 1);
	uint Tenth = uint(Conf*10);
	asserta(Tenth >= 0 && Tenth <= 10);
	static const char Symbols[12] = "___.,/:=@*^";
	//                               01234567890
	return Symbols[Tenth];
	}

static void Do1(FILE *fOut, const Ensemble &E, uint MSAIndex,
  const string &ConfLabel, int Dec)
	{
	const MSA &M = E.GetMSA(MSAIndex);
	const uint ColCount = M.GetColCount();
	string ConfSeq;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		double Conf = E.GetConf_MSACol(MSAIndex, ColIndex);
		char c = '?';
		switch (Dec)
			{
		case -1:
			c = ConfToChar_1(Conf);
			break;

		case 1:
			c = ConfToChar1(Conf);
			break;
		case 2:
			c = ConfToChar2(Conf);
			break;
		default:
			asserta(false);
			}
		ConfSeq += c;
		}

	Pf(fOut, ">%s\n", ConfLabel.c_str());
	Pf(fOut, "%s\n", ConfSeq.c_str());
	}

void cmd_addconfseq()
	{
	const string &InputFileName = opt(addconfseq);
	const string RefFileName = opt(ref);
	const string &OutputFileName = opt(output);

	string ConfLabel = "_conf_";
	if (optset_label)
		ConfLabel = opt(label);

	MSA Ref;
	if (optset_ref)
		Ref.FromFASTAFile_PreserveCase(RefFileName);

	Ensemble E;
	E.FromFile(InputFileName);
	if (optset_ref)
		E.SortMSA(Ref);

	FILE *fOut = CreateStdioFile(OutputFileName);

	const uint MSACount = E.GetMSACount();
	const uint SeqCount = E.GetSeqCount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &MSAName = E.GetMSAName(MSAIndex);
		Pf(fOut, "<%s\n", MSAName.c_str());
		const MSA &M = E.GetMSA(MSAIndex);

		const uint ColCount = M.GetColCount();
		const uint MSASeqCount = M.GetSeqCount();
		asserta(MSASeqCount == SeqCount);

		if (opt(confseq1))
			Do1(fOut, E, MSAIndex, ConfLabel, -1);
		else
			{
			Do1(fOut, E, MSAIndex, ConfLabel, 1);
			Do1(fOut, E, MSAIndex, ConfLabel + "2", 2);
			}
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const char *S = M.m_szSeqs[SeqIndex];
			const char *Label = M.m_szNames[SeqIndex];
			Pf(fOut, ">%s\n", Label);
			Pf(fOut, "%*.*s\n", ColCount, ColCount, S);
			}
		}

	CloseStdioFile(fOut);
	}
