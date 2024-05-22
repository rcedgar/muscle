#include "muscle.h"

/***
Input is "d:/a/res/newbench/out/selectpfams.tsv" from -newbench_selectpfams

      0           1       2     3        4       5       6
Uniprot          Sp       L     N       PF      Lo      Hi...
Q9IBF7  OTOMP_ONCMY     367     1       PF00405 27      363
P0C2M5  PSTS2_STRPN     291     2       PF12849 4       155     PF12849 172     290
Q26643  TRF_SARPE       629     3       PF00405 27      362     PF00405 372     474     PF00405 498     627

Output is sorted by DomStr e.g. PF00009+PF03144+PF03143, given sequence
  appears exactly once in the output (no duplicates).
All sequences with given DomStr are globally alignable.
***/

static map<string, vector<string> > g_PFToMultiDomStrs;
static map<string, vector<string> > g_DomStrToLines;
static map<string, uint> g_UpToClusterIndex;
static MultiSequence g_Seqs;
static map<string, uint> g_UpToSeqIndex;
static string g_OutDirFa = "d:/int/newbench/fa/";
static string g_OutDirTsv= "d:/int/newbench/tsv/";
static set<string> g_DonePFs;

static const uint MIN_SIZE = 10;

static bool ContainsPFExactlyOnce(const string &DomStr, const string &PF)
	{
	vector<string> Fields;
	Split(DomStr, Fields, '_');
	uint n = 0;
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		if (Fields[i] == PF)
			++n;
		}
	return n == 1;
	}

static bool DomStrHasDonePF(const string &DomStr)
	{
	vector<string> PFs;
	Split(DomStr, PFs, '_');
	for (uint i = 0; i < SIZE(PFs); ++i)
		if (g_DonePFs.find(PFs[i]) != g_DonePFs.end())
			return true;
	return false;
	}

static void OutputDomStr(FILE *fFa, FILE *fTsv, const string &DomStr, set<string> &DoneLabels)
	{
	asserta(g_DomStrToLines.find(DomStr) != g_DomStrToLines.end());
	const vector<string> &Lines = g_DomStrToLines[DomStr];
	const uint SeqCountThisDomStr = SIZE(Lines);
	for (uint SeqIndexThisDomStr = 0; SeqIndexThisDomStr < SeqCountThisDomStr; ++SeqIndexThisDomStr)
		{
		const string &Line = Lines[SeqIndexThisDomStr];
		string Up;
		uint n = SIZE(Line);
		for (uint i = 0; i < n; ++i)
			{
			char c = Line[i];
			if (c == '\t')
				break;
			Up += c;
			}

		if (DoneLabels.find(Up) != DoneLabels.end())
			continue;
		DoneLabels.insert(Up);

		map<string, uint>::const_iterator iter = g_UpToSeqIndex.find(Up);
		asserta(iter != g_UpToSeqIndex.end());
		uint SeqIndex = iter->second;
		const char *Seq = g_Seqs.GetCharPtr(SeqIndex);
		const uint L = g_Seqs.GetSeqLength(SeqIndex);
		if (fFa != 0)
			SeqToFasta(fFa, Seq, L, Up.c_str());
		if (fTsv != 0)
			fprintf(fTsv, "%s\t%s\n", DomStr.c_str(), Line.c_str());
		}
	}

static void DoDomStrs(const string &PF, const set<string> &ArgDomStrs)
	{
	if (g_DonePFs.find(PF) != g_DonePFs.end())
		return;

	vector<string> DomStrs;
	for (set<string>::const_iterator iter = ArgDomStrs.begin();
	  iter != ArgDomStrs.end(); ++iter)
		{
		const string &DomStr = *iter;
		if (DomStrHasDonePF(DomStr))
			continue;
		if (!ContainsPFExactlyOnce(DomStr, PF))
			continue;
		DomStrs.push_back(DomStr);
		}
	if (SIZE(DomStrs) < 2)
		return;

	uint Size = 0;
	for (uint i = 0; i < SIZE(DomStrs); ++i)
		{
		const string &DomStr = DomStrs[i];
		const vector<string> &Lines = g_DomStrToLines[DomStr];
		Size += SIZE(Lines);
		}
	if (Size < MIN_SIZE)
		return;

	const string Name = PF + "_local";
	string OutputFastaFileName = g_OutDirFa + Name;
	string OutputTsvFileName = g_OutDirTsv + Name;
	FILE *fFa = CreateStdioFile(OutputFastaFileName);
	FILE *fTsv = CreateStdioFile(OutputTsvFileName);

	Log("DoDomStrs(PF=%s)\n", PF.c_str());
	set<string> DoneLabels;
	for (uint i = 0; i < SIZE(DomStrs); ++i)
		{
		const string &DomStr = DomStrs[i];
		Log("  DomStr %s\n", DomStr.c_str());
		for (uint i = 0; i < SIZE(DomStrs); ++i)
			OutputDomStr(fFa, fTsv, DomStr, DoneLabels);
		}

	set<string> PFSet;
	for (uint i = 0; i < SIZE(DomStrs); ++i)
		{
		const string &DomStr = DomStrs[i];
		vector<string> PFs;
		Split(DomStr, PFs, '_');
		for (uint i = 0; i < SIZE(PFs); ++i)
			{
			const string &PF2 = PFs[i];
			PFSet.insert(PF2);
			}
		}

	for (set<string>::const_iterator p = PFSet.begin();
	  p != PFSet.end(); ++p)
		{
		const string &PF2 = *p;
		asserta(g_DonePFs.find(PF2) == g_DonePFs.end());
		g_DonePFs.insert(PF2);
		Log("  done %s\n", PF2.c_str());
		}

	CloseStdioFile(fFa);
	CloseStdioFile(fTsv);
	}

static void DoDomStr(const string &DomStr, uint MinSize)
	{
	asserta(g_DomStrToLines.find(DomStr) != g_DomStrToLines.end());
	const vector<string> &Lines = g_DomStrToLines[DomStr];
	const uint SeqCountThisDomStr = SIZE(Lines);

	const uint LineCount = SIZE(Lines);
	if (SeqCountThisDomStr < MinSize)
		return;

	string OutputFastaFileName = g_OutDirFa + DomStr;
	string OutputTsvFileName = g_OutDirTsv + DomStr;
	FILE *fFa = CreateStdioFile(OutputFastaFileName);
	FILE *fTsv = CreateStdioFile(OutputTsvFileName);
	set<string> DoneLabels;
	OutputDomStr(fFa, fTsv, DomStr, DoneLabels);
	CloseStdioFile(fFa);
	CloseStdioFile(fTsv);
	}

static void ParseLine(const string &Line, string &DomStr,
  vector<string> &PFs)
	{
	DomStr.clear();
	PFs.clear();
	vector<string> Fields;
	Split(Line, Fields, '\t');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount >= 7);
	const uint PFCount = StrToUint(Fields[3]);
	asserta(FieldCount == 3*PFCount + 4);
	for (uint i = 0; i < PFCount; ++i)
		{
		const string &PF = Fields[4+3*i];
		asserta(StartsWith(PF, "PF"));
		PFs.push_back(PF);
		if (i > 0)
			DomStr += "_";
		DomStr += PF;
		}
	}

void cmd_newbench_pfamgroups()
	{
	string Line;
	vector<string> Fields;

	g_Seqs.FromFASTA(opt(input));
	const uint SeqCount = g_Seqs.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Up = g_Seqs.GetLabelStr(i);
		g_UpToSeqIndex[Up] = i;
		}

	FILE *f = OpenStdioFile(g_Arg1); // "d:/a/res/newbench/out/selectpfams.tsv"
	vector<string> Lines;

	map<string, set<string> > PFToDomStrs;
	uint MaxDomStrsPerPF = 0;
	while (ReadLineStdioFile(f, Line))
		{
		Lines.push_back(Line);
		string DomStr;
		vector<string> PFs;
		ParseLine(Line, DomStr, PFs);
		uint n = SIZE(PFs);
		for (uint i = 0; i < n; ++i)
			{
			const string &PF = PFs[i];
			if (PFToDomStrs.find(PF) == PFToDomStrs.end())
				{
				set<string> Empty;
				PFToDomStrs[PF] = Empty;
				}
			PFToDomStrs[PF].insert(DomStr);
			const uint DomStrCount = SIZE(PFToDomStrs[PF]);
			MaxDomStrsPerPF = max(DomStrCount, MaxDomStrsPerPF);
			}

		if (g_DomStrToLines.find(DomStr) == g_DomStrToLines.end())
			{
			vector<string> Empty;
			g_DomStrToLines[DomStr] = Empty;
			}
		g_DomStrToLines[DomStr].push_back(Line);
		}

	for (uint DomStrsPerPF = MaxDomStrsPerPF; DomStrsPerPF > 1; --DomStrsPerPF)
		{
		for (map<string, set<string> >::const_iterator iter = PFToDomStrs.begin();
		  iter != PFToDomStrs.end(); ++iter)
			{
			const string &PF = iter->first;
			const set<string> &DomStrs = iter->second;
			if (SIZE(DomStrs) != DomStrsPerPF)
				continue;
			DoDomStrs(PF, DomStrs);
			}
		}

	for (map<string, vector<string> >::const_iterator iter = g_DomStrToLines.begin();
	  iter != g_DomStrToLines.end(); ++iter)
		{
		set<string> PFSet;
		const string &DomStr = iter->first;
		vector<string> PFs;
		Split(DomStr, PFs, '_');
		bool Done = false;
		for (uint i = 0; i < SIZE(PFs); ++i)
			{
			const string &PF2 = PFs[i];
			if (g_DonePFs.find(PF2) != g_DonePFs.end())
				{
				Done = true;
				break;
				}
			PFSet.insert(PF2);
			}
		if (Done)
			continue;

		DoDomStr(DomStr, MIN_SIZE);

		Log("DomStr(%s)\n", DomStr.c_str());
		for (set<string>::const_iterator p = PFSet.begin();
		  p != PFSet.end(); ++p)
			{
			const string &PF2 = *p;
			asserta(g_DonePFs.find(PF2) == g_DonePFs.end());
			g_DonePFs.insert(PF2);
			Log("  done %s\n", PF2.c_str());
			}
		}
	}
