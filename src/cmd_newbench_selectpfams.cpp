#include "muscle.h"
#include "sort.h"

bool GetMSAColIsAligned(const MSA &Aln, uint Col);
void GetMSAColAlignedVec(const MSA &Aln,
  vector<bool> &AlignedVec, ptr_GetMSAColIsAligned pFn);
MSA *SqueezeInserts(const MSA &Aln, ptr_GetMSAColIsAligned pFn);

static const char *PFAMRegionsTsvFN =
  "d:/int/afdb_swissprot_bali/out/pfam_regions_intersect.tsv";

static const uint MIN_SEQ_COUNT = 8;
static const uint MIN_CORE_COL_COUNT = 10;
static const double MIN_PRIMARY_CORE_FRACT = 0.666;
static const double MIN_DOMAIN_LENGTH = 30;
static const double MIN_COVERAGE_GLOBAL = 0.8;
static const double MIN_DOM_FRACT = 0.1;

static uint g_NOK = 0;
static FILE *g_fRep = 0;

/***
head /d/int/afdb_swissprot_bali/out/pfam_regions_intersect.tsv
0       1               2       3       4       5
Uniprot Swissprot       PF      Lo		Hi		Length
P25344  STE50_YEAST     PF09235 30      104     346
P15306  TRBM_MOUSE      PF09064 405     438     577
Q71U07  TRBM_SAISC      PF09064 406     439     575
Q5W7P8  TRBM_CANLF      PF09064 407     440     578
P41156  ETS1_RAT        PF02198 53      136     441
P27577  ETS1_MOUSE      PF02198 53      136     440
Q93Z38  TAR4_ARATH      PF04863 33      89      463
***/

// d:/int/afdb_swissprot_bali/pfams/PF00012/mustang.afa
// d:/int/afdb_swissprot_bali/pfams/PF00012/id90.fa
static string g_PfamsDir = "d:/int/afdb_swissprot_bali/pfams/";
static string g_OutDir = "d:/int/newbench/";

static string g_PF;
static uint g_PFix;

static string g_AnnotRow;

// one entry for each unique PF
static vector<string> g_UniquePFs;
static map<string, uint> g_PFToPFix;
static vector<vector<uint> > g_PFixToRegionixs;

// one entry for each unique Uniprot accession
static vector<string> g_UniqueUps;
static map<string, uint> g_UpToUpix;
static vector<vector<uint> > g_UpixToRegionixs;
static map<string, string> g_UpToSp;

// one entry for each region
static vector<string> g_Ups;
static vector<string> g_Sps;
static vector<string> g_PFs;
static vector<uint> g_Los;
static vector<uint> g_His;
static vector<uint> g_Ls;

// one entry per sequence for current Aln
static vector<string> g_AlnUps;
static vector<string> g_AlnSps;
static vector<vector<uint> > g_AlnRegionixVec;
static vector<vector<uint> > g_PosToColVec;
static vector<vector<uint> > g_PosToPFixVec;
static vector<vector<uint> > g_ColToPFixVec;
static vector<string> g_UngappedSeqs;
static vector<string> g_DomStrs;
static vector<double> g_DomCoverages;
static vector<string> g_BeforePrimaryPFs;
static vector<string> g_AfterPrimaryPFs;

// For current Aln
static vector<uint> g_AlnPFixs;
static map<uint, char> g_AlnPFixToChar;
static vector<string> g_UniqueDomStrs;
static vector<uint> g_UniqueDomStrCounts;
static map<string, uint> g_DomStrToCount;
static vector<double> g_UniqueDomMinCoverages;
static double g_Globalness;
static uint g_CoreCount;
static uint g_CorePrimaryCount;
static uint g_CoreOtherCount;

static uint GetUngappedSeqLength(uint SeqIndex)
	{
	asserta(SeqIndex < SIZE(g_UngappedSeqs));
	uint L = SIZE(g_UngappedSeqs[SeqIndex]);
	return L;
	}

static uint GetPFix(const string &PF, bool NewOk)
	{
	map<string, uint>::const_iterator iter =
	  g_PFToPFix.find(PF);
	uint PFix = UINT_MAX;
	if (iter == g_PFToPFix.end())
		{
		if (!NewOk)
			Die("GetPFix(%s) not found", PF.c_str());

		PFix = SIZE(g_UniquePFs);
		g_UniquePFs.push_back(PF);
		g_PFToPFix[PF] = PFix;

		vector<uint> Empty;
		asserta(SIZE(g_PFixToRegionixs) == PFix);
		g_PFixToRegionixs.push_back(Empty);
		}
	else
		PFix = iter->second;
	return PFix;
	}

static uint GetUpix(const string &Up, bool NewOk)
	{
	map<string, uint>::const_iterator iter =
	  g_UpToUpix.find(Up);
	uint Upix = UINT_MAX;
	if (iter == g_UpToUpix.end())
		{
		if (!NewOk)
			Die("GetUpix(%s) not found", Up.c_str());

		Upix = SIZE(g_UniqueUps);
		asserta(SIZE(g_UpixToRegionixs) == Upix);
		g_UniqueUps.push_back(Up);
		g_UpToUpix[Up] = Upix;

		vector<uint> Empty;
		g_UpixToRegionixs.push_back(Empty);
		}
	else
		Upix = iter->second;
	return Upix;
	}

static void ReadPFAMRegions()
	{
	FILE *f = OpenStdioFile(PFAMRegionsTsvFN);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 6);

		const string &Up = Fields[0];
		const string &Sp = Fields[1];
		const string &PF = Fields[2];
		uint Lo = StrToUint(Fields[3]);
		uint Hi = StrToUint(Fields[4]);
		uint L = StrToUint(Fields[5]);
		asserta(Lo < Hi);
		asserta(Hi <= L);

		map<string, string>::const_iterator iter = g_UpToSp.find(Up);
		if (iter == g_UpToSp.end())
			g_UpToSp[Up] = Sp;
		else
			asserta(iter->second == Sp);

		uint PFix = GetPFix(PF, true);
		uint Upix = GetUpix(Up, true);
		uint Regionix = SIZE(g_Ups);
		g_Ups.push_back(Up);
		g_Sps.push_back(Sp);
		g_PFs.push_back(PF);
		g_Los.push_back(Lo);
		g_His.push_back(Hi);
		g_Ls.push_back(L);

		asserta(Upix < SIZE(g_UpixToRegionixs));
		g_UpixToRegionixs[Upix].push_back(Regionix);
		}
	}

static bool HasPrimaryPFRepeat(const MSA &Aln)
	{
	const uint SeqCount = Aln.GetSeqCount();
	asserta(SIZE(g_AlnRegionixVec) == SeqCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint PrimaryCount = 0;
		const vector<uint> &Regionixs = g_AlnRegionixVec[SeqIndex];
		for (uint i = 0; i < SIZE(Regionixs); ++i)
			{
			uint Regionix = Regionixs[i];
			asserta(Regionix < SIZE(g_PFs));
			const string &PF = g_PFs[Regionix];
			if (PF == g_PF)
				++PrimaryCount;
			}
		asserta(PrimaryCount > 0);
		if (PrimaryCount > 1)
			{
			fprintf(g_fRep, "%s has primary repeat\n", Aln.GetLabel(SeqIndex));
			return true;
			}
		}
	return false;
	}

static double GetMeanPrimaryDomainLength(const MSA &Aln)
	{
	const uint SeqCount = Aln.GetSeqCount();
	asserta(SIZE(g_AlnRegionixVec) == SeqCount);
	uint SumPrimaryLength = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint PrimaryCount = 0;
		const vector<uint> &Regionixs = g_AlnRegionixVec[SeqIndex];
		for (uint i = 0; i < SIZE(Regionixs); ++i)
			{
			uint Regionix = Regionixs[i];
			asserta(Regionix < SIZE(g_PFs));
			const string &PF = g_PFs[Regionix];
			if (PF == g_PF)
				{
				uint Lo = g_Los[Regionix];
				uint Hi = g_His[Regionix];
				asserta(Lo < Hi);
				SumPrimaryLength += Hi - Lo + 1;
				}
			}
		}
	double Mean = SumPrimaryLength/SeqCount;
	return Mean;
	}

static void SetAlnRegions(const MSA &Aln)
	{
	g_AlnRegionixVec.clear();
	g_AlnPFixs.clear();
	g_AlnPFixToChar.clear();
	g_DomStrToCount.clear();
	g_DomStrs.clear();
	g_DomCoverages.clear();

	g_AlnPFixToChar[UINT_MAX] = '.';
	
	const uint SeqCount = Aln.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Up = (string) Aln.GetLabel(SeqIndex);
		map<string, uint>::const_iterator iter =
		  g_UpToUpix.find(Up);
		asserta(iter != g_UpToUpix.end());

		uint Upix = iter->second;
		const string &Up2 = g_UniqueUps[Upix];
		asserta(Up2 == Up);

		map<string, string>::const_iterator iter2 = g_UpToSp.find(Up);
		asserta(iter2 != g_UpToSp.end());
		const string &Sp = iter2->second;

		asserta(Upix < SIZE(g_UpixToRegionixs));
		const vector<uint> &Regionixs = g_UpixToRegionixs[Upix];
		g_AlnRegionixVec.push_back(Regionixs);

		string DomStr;
		const uint n = SIZE(Regionixs);
		fprintf(g_fRep, "Regions[%4u] %s %s %u:", SeqIndex, Up.c_str(), Sp.c_str(), n);
		uint SumDomLens = 0;
		for (uint i = 0; i < n; ++i)
			{
			uint Regionix = Regionixs[i];
			asserta(Regionix < SIZE(g_Ups));
			asserta(Regionix < SIZE(g_Sps));
			asserta(Regionix < SIZE(g_PFs));
			asserta(Regionix < SIZE(g_Los));
			asserta(Regionix < SIZE(g_His));
			asserta(Regionix < SIZE(g_Ls));

			const string &Up3 = g_Ups[Regionix];
			const string &Sp = g_Sps[Regionix];
			const string &PF = g_PFs[Regionix];
			uint Lo = g_Los[Regionix];
			uint Hi = g_His[Regionix];
			uint L = g_Ls[Regionix];
			uint DomLen = Hi - Lo + 1;
			if (double(DomLen)/L >= MIN_DOM_FRACT)
				{
				SumDomLens += DomLen;
				if (!DomStr.empty())
					DomStr += "+";
				DomStr += PF;
				}
			uint PFix = GetPFix(PF, false);

			if (g_AlnPFixToChar.find(PFix) == g_AlnPFixToChar.end())
				{
				uint i = SIZE(g_AlnPFixs);
				g_AlnPFixs.push_back(PFix);
				if (PFix == g_PFix)
					g_AlnPFixToChar[PFix] = '@';
				else
					g_AlnPFixToChar[PFix] = 'A' + i;
				}
			
			asserta(Up3 == Up);
			fprintf(g_fRep, " %u/%s:%u-%u(%u)", PFix, PF.c_str(), Lo, Hi, L);
			}
		fprintf(g_fRep, "  >%s\n", Aln.GetLabel(SeqIndex));

		if (DomStr == "")
			DomStr = "-";
		uint L = GetUngappedSeqLength(SeqIndex);
		if (SumDomLens >= L)
			SumDomLens = L;
		double DomCoverage = double(SumDomLens)/L;

		g_DomStrs.push_back(DomStr);
		g_DomCoverages.push_back(DomCoverage);
		if (g_DomStrToCount.find(DomStr) == g_DomStrToCount.end())
			g_DomStrToCount[DomStr] = 1;
		else
			g_DomStrToCount[DomStr] += 1;
		}

	asserta(SIZE(g_DomStrs) == SeqCount);
	asserta(SIZE(g_DomCoverages) == SeqCount);

	fprintf(g_fRep, "\n");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		fprintf(g_fRep, "DomStr[%4u]", SeqIndex);
		fprintf(g_fRep, "  %s (%.4f)\n",
		  g_DomStrs[SeqIndex].c_str(), g_DomCoverages[SeqIndex]);
		}

	VecToCountMap(g_DomStrs, g_DomStrToCount);
	CountMapToVecs(g_DomStrToCount, g_UniqueDomStrs, g_UniqueDomStrCounts);
	uint N = SIZE(g_UniqueDomStrs);
	vector<uint> Order(SeqCount);
	QuickSortOrderDesc(g_UniqueDomStrCounts.data(), N, Order.data());

	g_UniqueDomMinCoverages.clear();
	fprintf(g_fRep, "\n");
	bool GlobalFound = false;
	for (uint k = 0; k < N; ++k)
		{
		uint i = Order[k];
		const string &DomStr = g_UniqueDomStrs[i];
		uint Count = g_UniqueDomStrCounts[i];
		double MinCvg = DBL_MAX;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			if (g_DomStrs[SeqIndex] == DomStr)
				{
				double DomCvg = g_DomCoverages[SeqIndex];
				MinCvg = min(MinCvg, DomCvg);
				}
			}
		g_UniqueDomMinCoverages.push_back(MinCvg);
		fprintf(g_fRep, "DSN[%u]  %s  %6.4f  %u", k, DomStr.c_str(), MinCvg, Count);
		if (k == 0 && N == 1 && MinCvg >= MIN_COVERAGE_GLOBAL)
			{
			GlobalFound = true;
			fprintf(g_fRep, " << global");
			}
		fprintf(g_fRep, "\n");
		}

	uint Top = Order[0];
	asserta(Top < SIZE(g_UniqueDomStrs));
	asserta(Top < SIZE(g_UniqueDomStrCounts));
	asserta(Top < SIZE(g_UniqueDomMinCoverages));
	const string &TopDomStr = g_UniqueDomStrs[Top];
	uint TopCount = g_UniqueDomStrCounts[Top];
	asserta(TopCount <= SeqCount);
	double TopFract = double(TopCount)/SeqCount;
	double TopMinCvg = g_UniqueDomMinCoverages[Top];
	g_Globalness = TopFract*TopMinCvg;
	fprintf(g_fRep, "TopDS %s: count %u, fract %.4f, mincvg %.4f, GF %c, score %.3g\n",
	  TopDomStr.c_str(), TopCount, TopFract, TopMinCvg, yon(GlobalFound), g_Globalness);
	}

static void GetBeforeAfterPFs(const MSA &Aln, uint SeqIndex, 
  string &Before, string &After)
	{
	Before = ".";
	After = ".";
	asserta(SeqIndex < SIZE(g_AlnRegionixVec));
	const vector<uint> &Regionixs = g_AlnRegionixVec[SeqIndex];
	const uint n = SIZE(Regionixs);
	for (uint i = 0; i < n; ++i)
		{
		uint Regionix = Regionixs[i];
		asserta(Regionix < SIZE(g_PFs));
		const string &RegionPF = g_PFs[Regionix];
		if (RegionPF == g_PF)
			{
			if (i > 0)
				{
				uint BeforeRegionix = Regionixs[i-1];
				asserta(BeforeRegionix < SIZE(g_PFs));
				Before = g_PFs[BeforeRegionix];
				}
			if (i + 1 < n)
				{
				uint AfterRegionix = Regionixs[i+1];
				asserta(AfterRegionix < SIZE(g_PFs));
				After = g_PFs[AfterRegionix];
				}
			}
		}
	}

static void SetBeforeAfterPFs(const MSA &Aln)
	{
	g_BeforePrimaryPFs.clear();
	g_AfterPrimaryPFs.clear();
	const uint SeqCount = Aln.GetSeqCount();
	fprintf(g_fRep, "       %10.10s  %10.10s\n", "Before", "After");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Before, After;
		GetBeforeAfterPFs(Aln, SeqIndex, Before, After);
		g_BeforePrimaryPFs.push_back(Before);
		g_AfterPrimaryPFs.push_back(After);
	
		fprintf(g_fRep, "[%4u] %10.10s  %10.10s  >%s\n",
		  SeqIndex, Before.c_str(), After.c_str(), Aln.GetLabel(SeqIndex));
		}

	map<string, uint> CountMap;
	VecToCountMap(g_BeforePrimaryPFs, CountMap);
	vector<string> BeforePFs;
	vector<string> AfterPFs;
	vector<uint> BeforeCounts;
	vector<uint> AfterCounts;
	CountMapToVecs(CountMap, BeforePFs, BeforeCounts);
	CountMapToVecs(CountMap, AfterPFs, AfterCounts);
	uint BeforeTopCount = BeforeCounts[0];
	uint AfterTopCount = AfterCounts[0];
	double Agreement = 
	  double(BeforeTopCount)*double(AfterTopCount)/(SeqCount*SeqCount);
	fprintf(g_fRep, "BeforeAfterAgreement = %.4f\n", Agreement);
	}

static void SetAlnLabels(MSA &Aln)
	{
	g_AlnUps.clear();
	g_AlnSps.clear();

	const uint SeqCount = Aln.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Label = (string) Aln.GetLabel(SeqIndex);
		asserta(EndsWith(Label, ".pdb"));
		string Up = Label.substr(0, Label.size() - 4);
		map<string, uint>::const_iterator iter =
		  g_UpToUpix.find(Up);
		if (iter == g_UpToUpix.end())
			{
			fprintf(g_fRep, "Not found >%s\n", Up.c_str());
			continue;
			}
		Aln.m_szNames[SeqIndex] = mystrsave(Up.c_str());

		uint Upix = iter->second;
		const string &Up2 = g_UniqueUps[Upix];
		asserta(Up2 == Up);

		map<string, string>::const_iterator iter2 = g_UpToSp.find(Up);
		asserta(iter2 != g_UpToSp.end());
		const string &Sp = iter2->second;

		g_AlnUps.push_back(Up);
		g_AlnSps.push_back(Sp);
		}
	}

// Upper case if all letters, otherwise lower
static void SetAlnCase(MSA &Aln)
	{
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	asserta(SIZE(g_AnnotRow) == ColCount);
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		bool HasGap = Aln.ColumnHasGap(ColIndex);
		char c = g_AnnotRow[ColIndex];
		bool IsCoreCol = (c != '~');
		if (IsCoreCol)
			asserta(!HasGap);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char c = Aln.GetChar(SeqIndex, ColIndex);
			if (isalpha(c))
				c = (IsCoreCol ? toupper(c) : tolower(c));
			Aln.SetChar(SeqIndex, ColIndex, c);
			}
		}
	}

static void SetPosToColVec(const MSA &Aln)
	{
	g_PosToColVec.clear();
	const uint SeqCount = Aln.GetSeqCount();
	g_PosToColVec.resize(SeqCount);
	for (uint i = 0; i < SeqCount ; ++i)
		Aln.GetPosToCol(i, g_PosToColVec[i]);
	}

static void SetPosToPFixVec(const MSA &Aln)
	{
	g_PosToPFixVec.clear();
	const uint SeqCount = Aln.GetSeqCount();
	g_PosToPFixVec.resize(SeqCount);
	uint OverlapCount = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(g_AlnRegionixVec));
		uint L = GetUngappedSeqLength(SeqIndex);
		vector<uint> &PosToPFix = g_PosToPFixVec[SeqIndex];
		PosToPFix.resize(L, UINT_MAX);

		const vector<uint> &Regionixs = g_AlnRegionixVec[SeqIndex];
		const uint n = SIZE(Regionixs);
		for (uint i = 0; i < n; ++i)
			{
			uint Regionix = Regionixs[i];
			asserta(Regionix < SIZE(g_Ups));
			asserta(Regionix < SIZE(g_Sps));
			asserta(Regionix < SIZE(g_PFs));
			asserta(Regionix < SIZE(g_Los));
			asserta(Regionix < SIZE(g_His));
			asserta(Regionix < SIZE(g_Ls));

			const string &Up3 = g_Ups[Regionix];
			const string &Sp = g_Sps[Regionix];
			const string &PF = g_PFs[Regionix];
			uint Lo = g_Los[Regionix];
			uint Hi = g_His[Regionix];
			uint L2 = g_Ls[Regionix];
			asserta(L2 == L);
			uint PFix = GetPFix(PF, false);
			asserta(Lo > 0);
			asserta(Lo < Hi);
			asserta(Hi <= L);
			for (uint Pos = Lo; Pos <= Hi; ++Pos)
				{
				asserta(Pos > 0 && Pos <= L);
				uint PFix2 = PosToPFix[Pos-1];
				if (PFix2 != UINT_MAX)
					{
					++OverlapCount;
					if (PFix2 != g_PFix)
						PosToPFix[Pos-1] = PFix;
					}
				else
					PosToPFix[Pos-1] = PFix;
				}
			}
		}
	if (OverlapCount > 0)
		fprintf(g_fRep, "%u region overlaps", OverlapCount);
	}

static void SetColToPFixVec(const MSA &Aln)
	{
	g_ColToPFixVec.clear();
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	g_ColToPFixVec.resize(SeqCount);
	fprintf(g_fRep, "\n");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		asserta(SeqIndex < SIZE(g_AlnRegionixVec));
		uint L = GetUngappedSeqLength(SeqIndex);
		const vector<uint> &PosToPFix = g_PosToPFixVec[SeqIndex];
		const vector<uint> &PosToCol = g_PosToColVec[SeqIndex];
		asserta(SIZE(PosToPFix) == L);
		asserta(SIZE(PosToCol) == L);
		vector<uint> &ColToPFix = g_ColToPFixVec[SeqIndex];
		ColToPFix.resize(ColCount, UINT_MAX);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Col = PosToCol[Pos];
			asserta(Col < ColCount);
			uint PFix = PosToPFix[Pos];
			ColToPFix[Col] = PFix;
			}

		for (uint Col = 0; Col < ColCount; ++Col)
			{
			uint PFix = ColToPFix[Col];
			map<uint, char>::const_iterator iter = g_AlnPFixToChar.find(PFix);
			asserta(iter != g_AlnPFixToChar.end());
			char c = iter->second;
			fprintf(g_fRep, "%c", c);
			}
		fprintf(g_fRep, " [%4u] ", SeqIndex);
		fprintf(g_fRep, " >%s\n", Aln.GetLabel(SeqIndex));
		}

	if (0)
		{
		fprintf(g_fRep, "\n");
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			asserta(SeqIndex < SIZE(g_AlnRegionixVec));
			uint L = GetUngappedSeqLength(SeqIndex);
			const vector<uint> &PosToPFix = g_PosToPFixVec[SeqIndex];
			asserta(SIZE(PosToPFix) == L);

			for (uint Pos = 0; Pos < L; ++Pos)
				{
				uint PFix = PosToPFix[Pos];
				map<uint, char>::const_iterator iter = g_AlnPFixToChar.find(PFix);
				asserta(iter != g_AlnPFixToChar.end());
				char c = iter->second;
				fprintf(g_fRep, "%c", c);
				}
			fprintf(g_fRep, " [%u] ", SeqIndex);
			fprintf(g_fRep, " >%s\n", Aln.GetLabel(SeqIndex));
			}
		}
	}

static void SetUngappedSeqs(const MSA &Aln)
	{
	g_UngappedSeqs.clear();
	const uint SeqCount = Aln.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Seq;
		Aln.GetUngappedSeqStr(SeqIndex, Seq);
		g_UngappedSeqs.push_back(Seq);
		}
	}

static uint GetConsensusPF(const MSA &Aln, uint Col)
	{
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	asserta(Col < ColCount);
	asserta(SIZE(g_ColToPFixVec) == SeqCount);
	uint ConsensusPFix = UINT_MAX;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const vector<uint> &ColToPFix = g_ColToPFixVec[SeqIndex];
		if (SIZE(ColToPFix) != ColCount)
			Die("SIZE(ColToPFix[%u])=%d, ColCount=%d",
			  SeqIndex, SIZE(ColToPFix), ColCount);
		uint PFix = ColToPFix[Col];
		if (PFix == UINT_MAX)
			return UINT_MAX;
		if (SeqIndex == 0)
			ConsensusPFix = PFix;
		else if (PFix != ConsensusPFix)
			return UINT_MAX;
		}
	return ConsensusPFix;
	}

// ^   no gaps, all primary PF
// :   no gaps. all non-primary PF
// ~   otherwise
static void SetColAnnots(const MSA &Aln)
	{
	g_CoreCount = 0;
	g_CorePrimaryCount = 0;
	g_CoreOtherCount = 0;

	g_AnnotRow.clear();
	const uint ColCount = Aln.GetColCount();
	uint CoreCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint PFix = GetConsensusPF(Aln, Col);
		if (PFix == UINT_MAX)
			{
			g_AnnotRow += '~';
			continue;
			}

		++CoreCount;
		if (PFix == g_PFix)
			{
			++g_CorePrimaryCount;
			g_AnnotRow +=  '^';
			}
		else
			{
			++g_CoreOtherCount;
			g_AnnotRow +=  ':';
			}
		}
	g_CoreCount = g_CorePrimaryCount + g_CoreOtherCount;
	}

static void GetPrimaryColLoHi(const MSA &Aln, uint SeqIndex, 
  uint &PrimaryColLo, uint &PrimaryColHi)
	{
	PrimaryColLo = UINT_MAX;
	PrimaryColHi = UINT_MAX;
	asserta(SeqIndex < SIZE(g_ColToPFixVec));
	const vector<uint> &ColToPFix = g_ColToPFixVec[SeqIndex];
	const uint ColCount = Aln.GetColCount();
	asserta(SIZE(ColToPFix) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint PFix = ColToPFix[Col];
		if (PFix == g_PFix)
			{
			if (PrimaryColLo == UINT_MAX)
				PrimaryColLo = Col;
			PrimaryColHi = Col;
			}
		}
	asserta(PrimaryColLo != UINT_MAX);
	asserta(PrimaryColHi != UINT_MAX);
	asserta(PrimaryColLo < PrimaryColHi);
	}

static MSA *Trim(const MSA &Aln)
	{
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	MSA *TmpAln = new MSA;
	TmpAln->Copy(Aln);
	asserta(TmpAln->GetSeqCount() == SeqCount);
	asserta(TmpAln->GetColCount() == ColCount);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint PrimaryColLo = UINT_MAX;
		uint PrimaryColHi = UINT_MAX;
		GetPrimaryColLoHi(Aln, SeqIndex, PrimaryColLo, PrimaryColHi);
		for (uint Col = 0; Col < PrimaryColLo; ++Col)
			TmpAln->SetChar(SeqIndex, Col, '.');
		for (uint Col = PrimaryColHi + 1 ; Col < ColCount; ++Col)
			TmpAln->SetChar(SeqIndex, Col, '.');
		}
	MSA *TrimmedAln = SqueezeInserts(*TmpAln, GetMSAColIsAligned);
	return TrimmedAln;
	}

//static void Do1(const string &SetName)
//	{
//	asserta(StartsWith(SetName, "PF"));
//	g_PF = SetName;
//	g_PFix = GetPFix(g_PF, true);
//
//	string MustangFN = g_PfamsDir + SetName + "/mustang.afa";
//	if (!StdioFileExists(MustangFN))
//		{
//		Warning("Not found %s", MustangFN.c_str());
//		return;
//		}
//	string ReportFN = g_OutDir + "report/" + SetName;
//	CloseStdioFile(g_fRep);
//	g_fRep = CreateStdioFile(ReportFN);
//
//	MSA Aln;
//	Aln.FromFASTAFile(MustangFN);
//	const uint SeqCount = Aln.GetSeqCount();
//	const uint ColCount = Aln.GetColCount();
//
//	if (SeqCount < MIN_SEQ_COUNT)
//		{
//		{
//		Log("%s rejected, %u seqs < %u\n",
//		  g_PF.c_str(), SeqCount, MIN_SEQ_COUNT);
//		fprintf(g_fRep, "%s rejected, %u seqs < %u\n",
//		  g_PF.c_str(), SeqCount, MIN_SEQ_COUNT);
//		return;
//		}
//		}
//
//	SetUngappedSeqs(Aln);
//	SetPosToColVec(Aln);
//	SetAlnLabels(Aln);
//	SetAlnRegions(Aln);
//	SetBeforeAfterPFs(Aln);
//	bool HasRepeat = HasPrimaryPFRepeat(Aln);
//	if (HasRepeat)
//		{
//		Log("%s rejected, primary repeat\n", g_PF.c_str());
//		fprintf(g_fRep, "rejected, primary repeat\n");
//		return;
//		}
//	double MeanLen = GetMeanPrimaryDomainLength(Aln);
//	if (MeanLen < MIN_DOMAIN_LENGTH)
//		{
//		Log("%s rejected, mean length %.1f < %.1f\n",
//		  g_PF.c_str(), MeanLen, MIN_DOMAIN_LENGTH);
//		fprintf(g_fRep, "rejected, mean length %.1f < %.1f\n",
//		  MeanLen, MIN_DOMAIN_LENGTH);
//		return;
//		}
//
//	SetPosToPFixVec(Aln);
//	SetColToPFixVec(Aln);
//	SetColAnnots(Aln);
//	if (g_CoreCount < MIN_CORE_COL_COUNT)
//		{
//		Log("%s rejected, not enough core cols (%u)\n", g_PF.c_str(), g_CoreCount);
//		fprintf(g_fRep, "rejected, not enough core cols (%u)\n", g_CoreCount);
//		return;
//		}
//	double PrimaryFract = double(g_CorePrimaryCount)/g_CoreCount;
//	if (PrimaryFract < MIN_PRIMARY_CORE_FRACT)
//		{
//		Log("%s rejected, primary core cols %u / %u = %.3g < %.3g\n",
//		  g_PF.c_str(), g_CorePrimaryCount, g_CoreCount, PrimaryFract, MIN_PRIMARY_CORE_FRACT);
//		fprintf(g_fRep, "%s rejected, primary core cols %u / %u = %.3g < %.3g\n",
//		  g_PF.c_str(), g_CorePrimaryCount, g_CoreCount, PrimaryFract, MIN_PRIMARY_CORE_FRACT);
//		return;
//		}
//
//	fprintf(g_fRep, "%s\n", g_AnnotRow.c_str());
//	for (map<uint, char>::const_iterator iter = g_AlnPFixToChar.begin();
//	  iter != g_AlnPFixToChar.end(); ++iter)
//		{
//		uint PFix = iter->first;
//		if (PFix == UINT_MAX || PFix == g_PFix)
//			continue;
//		char c = iter->second;
//		asserta(PFix < SIZE(g_UniquePFs));
//		const string &PF = g_UniquePFs[PFix];
//		fprintf(g_fRep, " %c=%s", c, PF.c_str());
//		}
//	fprintf(g_fRep, "\n");
//	fprintf(g_fRep, "%u / %u core cols\n", g_CoreCount, ColCount);
//
//	fprintf(g_fRep, "Mean primary domain length %.1f\n", MeanLen);
//
//	SetAlnCase(Aln);
//
//	MSA *SqueezedAln = SqueezeInserts(Aln, GetMSAColIsAligned);
//	string FN = g_OutDir + "ref/" + g_PF;
//	SqueezedAln->ToFASTAFile(FN);
//
//	MSA *TrimmedAln = Trim(Aln);
//	FN = g_OutDir + "ref_trimmed/" + g_PF;
//	TrimmedAln->ToFASTAFile(FN);
//
//	Log("%s OK\n", g_PF.c_str());
//	fprintf(g_fRep, "%s AcceptedOK\n", g_PF.c_str());
//	++g_NOK;
//	}
//
//void cmd_newbench_wrangle()
//	{
//	const string &SetsFN = g_Arg1;
//
//	ReadPFAMRegions();
//	
//	vector<string> SetNames;
//	ReadStringsFromFile(SetsFN, SetNames);
//
//	const uint SetCount = SIZE(SetNames);
//	for (uint SetIndex = 0; SetIndex < SetCount; ++SetIndex)
//		{
//		const string &SetName = SetNames[SetIndex];
//		ProgressStep(SetIndex, SetCount, "%s (%u OK)",
//		  g_PF.c_str(), g_NOK);
//		Do1(SetName);
//		}
//	}

static double GetAnnotCoverage(const vector<uint> &UpRegionIxs)
	{
	const uint n = SIZE(UpRegionIxs);
	asserta(n > 0);
	uint L = UINT_MAX;
	uint Upix = UINT_MAX;
	uint SumDomLens = 0;
	for (uint i = 0; i < n; ++i)
		{
		uint Regionix = UpRegionIxs[i];
		asserta(Regionix < SIZE(g_Los));
		asserta(Regionix < SIZE(g_His));
		asserta(Regionix < SIZE(g_Ls));

		const string &Up = g_Ups[Regionix];
		uint Upix2 = GetUpix(Up, false);
		const string &Sp = g_Sps[Regionix];
		const string &PF = g_PFs[Regionix];
		uint Lo = g_Los[Regionix];
		uint Hi = g_His[Regionix];
		uint L2 = g_Ls[Regionix];
		uint DomLen = Hi - Lo + 1;
		if (i == 0)
			{
			L = L2;
			Upix = Upix2;
			}
		else
			{
			asserta(L2 == L);
			asserta(Upix2 == Upix);
			}
		}
	return 0.0;
	}

static bool SelectUp(FILE *f, uint Upix)
	{
	const double MINCVG = 0.9;
	const double MINFRACT = 0.1;
	const string &Up = g_UniqueUps[Upix];
	asserta(Upix < SIZE(g_UpixToRegionixs));
	const vector<uint> &Regionixs = g_UpixToRegionixs[Upix];
	const uint N = SIZE(Regionixs);
	vector<string> PFs;
	vector<uint> Los;
	vector<uint> His;
	uint L = UINT_MAX;
	string Sp;
	uint SumDomLen = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint Regionix = Regionixs[i];
		const string &PF = g_PFs[Regionix];
		const string &Up2 = g_Ups[Regionix];
		const string &Sp2 = g_Sps[Regionix];
		uint Lo = g_Los[Regionix];
		uint Hi = g_His[Regionix];
		uint L2 = g_Ls[Regionix];

		if (i == 0)
			L = L2;
		else
			asserta(L2 == L);
		asserta(Up2 == Up);
		if (i == 0)
			Sp = Sp2;
		else
			assert(Sp2 == Sp);

	// 0-based in PFAM annotations
		asserta(Lo > 0);
		asserta(Lo < Hi);
		asserta(Hi <= L);
		uint DomLen = Hi - Lo + 1;
		double Fract = double(DomLen)/L;
		if (Fract >= MINFRACT)
			{
			for (uint j = 0; j < SIZE(PFs); ++j)
				{
				if (GetOverlap(Los[j], His[j], Lo, Hi) > 0)
					return false;
				}
			SumDomLen += DomLen;
			PFs.push_back(PF);
			Los.push_back(Lo);
			His.push_back(Hi);
			}
		}
	double Coverage = double(SumDomLen)/L;
	asserta(Coverage <= 1);
	if (Coverage < MINCVG)
		return false;
	if (f == 0)
		return true;

	const uint K = SIZE(Los);
	vector<uint> Order(K);
	QuickSortOrder(Los.data(), K, Order.data());

	fprintf(f, "%s", Up.c_str());
	fprintf(f, "\t%s", Sp.c_str());
	fprintf(f, "\t%u", L);
	fprintf(f, "\t%u", K);
	for (uint k = 0; k < K; ++k)
		{
		uint i = Order[k];
		uint Lo = Los[i];
		uint Hi = His[i];
		const string &PF = PFs[i];
		fprintf(f, "\t%s\t%u\t%u", PF.c_str(), Lo, Hi);
		}
	fprintf(f, "\n");
	return true;
	}

void cmd_newbench_selectpfams()
	{
	FILE *f = CreateStdioFile(g_Arg1);
	ReadPFAMRegions();
	const uint UpCount = SIZE(g_UniqueUps);
	uint SelectedCount = 0;
	for (uint Upix = 0; Upix < UpCount; ++Upix)
		{
		bool Selected = SelectUp(f, Upix);
		if (Selected)
			++SelectedCount;
		ProgressStep(Upix, UpCount, "%u selected (%.1f%%)", 
		  SelectedCount, GetPct(SelectedCount, Upix+1));
		}
	CloseStdioFile(f);

	ProgressLog("%u / %u selected (%.1f%%)", 
	  SelectedCount, UpCount, GetPct(SelectedCount, UpCount));
	}
