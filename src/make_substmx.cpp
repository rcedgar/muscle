#include "muscle.h"

static map<pair<string, string>, pair<string, string> > g_LabelPairToRows;

static void DeleteNotUpper(const string &RowIn, string &RowOut)
	{
	RowOut.clear();
	for (uint i = 0; i < SIZE(RowIn); ++i)
		{
		char c = RowIn[i];
		if (isupper(c) || c == '-')
			RowOut += c;
		}
	}

static uint g_AlnCount;
static vector<uint> g_PctIdCounts;
static vector<uint> g_LetterCounts;
static uint g_TotalLetters;
static vector<vector<uint> > g_LetterPairCounts;
static uint g_TotalPairs;
static uint g_MinPctId = 0;
static uint g_MaxPctId = 100;

static void AddPair(char a, char b)
	{
	if (a == '-' || b == 'b')
		return;
	uint Lettera = g_CharToLetterAmino[a];
	uint Letterb = g_CharToLetterAmino[b];
	if (Lettera >= 20 || Letterb >= 20)
		return;

	g_TotalLetters += 2;
	g_LetterCounts[Lettera] += 1;
	g_LetterCounts[Letterb] += 1;

	g_TotalPairs += 2;
	g_LetterPairCounts[Lettera][Letterb] += 1;
	g_LetterPairCounts[Letterb][Lettera] += 1;
	}

static void AddRows(const string &A, const string &B)
	{
	const uint ColCount = SIZE(A);
	asserta(SIZE(B) == ColCount);
	uint n = 0;
	uint N = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char a = A[Col];
		char b = B[Col];
		AddPair(a, b);
		if (isupper(a) || isupper(b))
			{
			++N;
			if (a == b)
				++n;
			}
		}
	if (N == 0)
		return;

	g_AlnCount += 1;
	uint PctId = (n*100)/N;
	asserta(PctId >= 0 && PctId <= 100);
	if (PctId < g_MinPctId || PctId > g_MaxPctId)
		return;
	g_PctIdCounts[PctId] += 1;

	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char a = A[Col];
		char b = B[Col];
		AddPair(a, b);
		}
	}

void cmd_make_substmx()
	{
	const string &InputFileName = opt(make_substmx);
	const string &OutputFileName = opt(output);
	string MxName = "MX";
	if (optset_label)
		MxName = opt(label);
	if (optset_minpctid)
		{
		asserta(optset_maxpctid);
		g_MinPctId = opt(minpctid);
		g_MaxPctId = opt(maxpctid);
		asserta(g_MinPctId <= g_MaxPctId);
		}

	vector<string> MSAFileNames;
	ReadStringsFromFile(InputFileName, MSAFileNames);

	const uint MSACount = SIZE(MSAFileNames);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &FileName = MSAFileNames[MSAIndex];
		ProgressStep(MSAIndex, MSACount, "Reading MSA %s", FileName.c_str());
		MultiSequence MSA;
		MSA.LoadMFA(FileName);
		const uint SeqCount = MSA.GetSeqCount();
		uint ColCount = MSA.GetColCount();
		for (uint i = 0; i < SeqCount; ++i)
			{
			const Sequence *Seqi = MSA.GetSequence(i);
			const string &Labeli = Seqi->GetLabel();
			string Rowi_;
			Seqi->GetSeqAsString(Rowi_);
			string Rowi;
			DeleteNotUpper(Rowi_, Rowi);
			for (uint j = 0; j < i; ++j)
				{
				const Sequence *Seqj = MSA.GetSequence(j);
				const string &Labelj = Seqj->GetLabel();
				string Rowj_;
				Seqj->GetSeqAsString(Rowj_);
				string Rowj;
				DeleteNotUpper(Rowj_, Rowj);
				asserta(SIZE(Rowi) == SIZE(Rowj));

				pair<string, string> Rows(Rowi, Rowj);
				pair<string, string> LabelPair(Labeli, Labelj);
				if (g_LabelPairToRows.find(LabelPair) == g_LabelPairToRows.end())
					g_LabelPairToRows[LabelPair] = Rows;
				else
					{
					const pair<string, string> &CurrRows = g_LabelPairToRows[LabelPair];
					if (CurrRows.first.size() < ColCount)
						g_LabelPairToRows[LabelPair] = Rows;
					}
				}
			}
		}

	g_PctIdCounts.resize(101);
	g_LetterCounts.resize(20, 0);
	g_LetterPairCounts.resize(20);
	for (uint i = 0; i < 20; ++i)
		g_LetterPairCounts[i].resize(20, 0);

	const uint AlnCount = SIZE(g_LabelPairToRows);
	uint Counter = 0;
	for (map<pair<string, string>, pair<string, string> >::const_iterator iter =
	  g_LabelPairToRows.begin(); iter != g_LabelPairToRows.end(); ++iter)
		{
		ProgressStep(Counter++, AlnCount, "Counting letters");
		const pair<string, string> &LabelPair = iter->first;
		const pair<string, string> &Rows = iter->second;
		//Log("\n");
		//Log("%s  %s\n", Rows.first.c_str(), LabelPair.first.c_str());
		//Log("%s  %s\n", Rows.second.c_str(), LabelPair.second.c_str());
		AddRows(Rows.first, Rows.second);
		}

	double Sum = 0;
	vector<double> Freqs(20);
	for (uint i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		uint n = g_LetterCounts[i];
		double Freq = double(n)/g_TotalLetters;
		Freqs[i] = Freq;
		Sum += Freq;
		Log("%c  %8.6f  %u\n", c, Freq, n);
		}
	Log("Sum = %8.6f\n", Sum);
	asserta(feq(Sum, 1.0));
	Log("\n\n");

	for (uint i = 0; i < 20; ++i)
		{
		Log("%c:  ", g_LetterToCharAmino[i]);
		for (uint j = 0; j < 20; ++j)
			Log("  %c=%10u", g_LetterToCharAmino[j], g_LetterPairCounts[i][j]);
		Log("\n");
		}

	Log("\n");
	Log("PctId distribution\n");
	for (int PctId = 100; PctId >= 0; --PctId)
		{
		uint n = g_PctIdCounts[PctId];
		double Freq = double(n)/g_AlnCount;
		Log("%d\t%u\t%.4g\n", PctId, n, Freq);
		}

	vector<vector<float> > PairFreqMx(20);
	for (uint i = 0; i < 20; ++i)
		PairFreqMx[i].resize(20, FLT_MAX);

	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			{
			uint nij = g_LetterPairCounts[i][j];
			double Freqij = double(nij)/g_TotalPairs;
			PairFreqMx[i][j] = (float) Freqij;
			PairFreqMx[j][i] = (float) Freqij;
			}
		}

	Log("\nPair frequencies\n ");
	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		Log("%10c", ci);
		}
	Log("\n");
	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		Log("%c", ci);
		for (uint j = 0; j <= i; ++j)
			{
			double Freqij = PairFreqMx[i][j];
			Log("%10.5f", Freqij);
			}
		Log("\n");
		}

	Log("\n");
	Log("Score matrix\n");
	Log(" ");
	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		Log("%10c", ci);
		}
	Log("\n");

	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		double Freqi = Freqs[i];
		Log("%c", ci);
		for (uint j = 0; j <= i; ++j)
			{
			double Freqj = Freqs[j];
			uint nij = g_LetterPairCounts[i][j];
			double Freqij = double(nij)/g_TotalPairs;
			double Score = log2(Freqij/(Freqi*Freqj));
			Log("%10.6f", Score);
			}
		Log("\n");
		}

	FILE *fOut = CreateStdioFile(OutputFileName);
	fprintf(fOut, "%s", MxName.c_str());
	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		fprintf(fOut, "\t%c", ci);
		}
	fprintf(fOut, "\n");

	for (uint i = 0; i < 20; ++i)
		{
		char ci = g_LetterToCharAmino[i];
		double Freqi = Freqs[i];
		fprintf(fOut, "%c", ci);
		for (uint j = 0; j < 20; ++j)
			{
			double Freqj = Freqs[j];
			uint nij = g_LetterPairCounts[i][j];
			double Freqij = double(nij)/g_TotalPairs;
			double Score = log2(Freqij/(Freqi*Freqj));
			fprintf(fOut, "\t%.3f", Score);
			}
		fprintf(fOut, "\n");
		}
	CloseStdioFile(fOut);
	}
