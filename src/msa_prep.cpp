#include "muscle.h"
#include "msa.h"
#include "msa_prep.h"
#include "multisequence.h"
#include "heatmapcolors.h"
#include <algorithm>
#include <cctype>

static bool MSAHasMixedCaseCols(const MSA &Aln)
	{
	const uint ColCount = Aln.GetColCount();
	const uint SeqCount = Aln.GetSeqCount();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint UpperCount = 0;
		uint LowerCount = 0;
		for (uint i = 0; i < SeqCount; ++i)
			{
			char c = Aln.GetChar(i, Col);
			if (isgap(c))
				continue;
			if (isupper(c))
				++UpperCount;
			else if (islower(c))
				++LowerCount;
			else
				Die("Unexpected sequence char '%c'", c);
			}
		if (UpperCount > 0 && LowerCount > 0)
			return true;
		}
	return false;
	}

static bool IsCoreCol_Occupancy(const MSA &Aln, uint Col, double MaxGapFract)
	{
	const uint SeqCount = Aln.GetSeqCount();
	uint GapCount = Aln.GetGapCount(Col);
	double GapFract = double(GapCount)/double(SeqCount);
	return GapFract <= MaxGapFract;
	}

static bool IsCoreCol_LowerCase(const MSA &Aln, uint Col)
	{
	const uint SeqCount = Aln.GetSeqCount();
	uint UpperCount = 0;
	uint LowerCount = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		char c = Aln.GetChar(i, Col);
		if (isgap(c))
			continue;
		if (isupper(c))
			++UpperCount;
		else if (islower(c))
			++LowerCount;
		else
			Die("Unexpected sequence char '%c'", c);
		}
	if (UpperCount > 0 && LowerCount > 0)
		Die("Mixed-case col %u", Col);
	return UpperCount > 0;
	}

static void GetCoreColMask(const MSA &Aln, const MSAPrepParams &Params,
  bool UseLowerCase, vector<bool> &IsCore)
	{
	IsCore.clear();
	const uint ColCount = Aln.GetColCount();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		bool Core = false;
		if (UseLowerCase)
			Core = IsCoreCol_LowerCase(Aln, Col);
		else
			Core = IsCoreCol_Occupancy(Aln, Col, Params.MaxGapFract);
		IsCore.push_back(Core);
		}
	}

static void EndTrimRange(const vector<bool> &IsCore, uint &Lo, uint &Hi)
	{
	const uint ColCount = SIZE(IsCore);
	Lo = 0;
	Hi = (ColCount == 0 ? 0 : ColCount - 1);
	while (Lo < ColCount && !IsCore[Lo])
		++Lo;
	while (Hi > Lo && !IsCore[Hi])
		--Hi;
	if (Lo >= ColCount || !IsCore[Lo])
		{
		Lo = 0;
		Hi = 0;
		if (ColCount == 0)
			Hi = UINT_MAX;
		}
	}

static void FindCoreBlocks(const vector<bool> &IsCore, uint TrimLo, uint TrimHi,
  uint MinLen, vector<pair<uint, uint> > &Blocks)
	{
	Blocks.clear();
	if (TrimHi < TrimLo)
		return;
	uint Start = UINT_MAX;
	for (uint Col = TrimLo; Col <= TrimHi; ++Col)
		{
		if (IsCore[Col])
			{
			if (Start == UINT_MAX)
				Start = Col;
			}
		else
			{
			if (Start != UINT_MAX)
				{
				if (Col - Start >= MinLen)
					Blocks.push_back(make_pair(Start, Col - 1));
				Start = UINT_MAX;
				}
			}
		}
	if (Start != UINT_MAX && TrimHi + 1 - Start >= MinLen)
		Blocks.push_back(make_pair(Start, TrimHi));
	}

struct InsertStats
	{
	uint NPresent = 0;
	uint MaxLen = 0;
	double FracPresent = 0;
	};

static InsertStats GetInsertStats(const MSA &Aln, uint Lo, uint Hi)
	{
	InsertStats Stats;
	const uint SeqCount = Aln.GetSeqCount();
	if (Lo > Hi)
		return Stats;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint Len = 0;
		for (uint Col = Lo; Col <= Hi; ++Col)
			{
			char c = Aln.GetChar(SeqIndex, Col);
			if (!isgap(c))
				++Len;
			}
		if (Len > 0)
			{
			++Stats.NPresent;
			if (Len > Stats.MaxLen)
				Stats.MaxLen = Len;
			}
		}
	if (SeqCount > 0)
		Stats.FracPresent = double(Stats.NPresent)/double(SeqCount);
	return Stats;
	}

static void MergeCoreBlocks(const MSA &Aln, vector<pair<uint, uint> > &Blocks,
  const MSAPrepParams &Params)
	{
	if (SIZE(Blocks) <= 1)
		return;
	vector<pair<uint, uint> > Merged;
	Merged.push_back(Blocks[0]);
	for (uint i = 1; i < SIZE(Blocks); ++i)
		{
		pair<uint, uint> &Prev = Merged.back();
		uint GapLo = Prev.second + 1;
		uint GapHi = Blocks[i].first - 1;
		bool DoMerge = false;
		if (GapHi >= GapLo)
			{
			uint GapWidth = GapHi - GapLo + 1;
			if (GapWidth <= Params.MergeGapCols)
				{
				InsertStats Stats = GetInsertStats(Aln, GapLo, GapHi);
				if (Stats.FracPresent <= Params.InsertMinorityFrac &&
				    Stats.MaxLen <= Params.InsertCapLen)
					DoMerge = true;
				}
			}
		else
			DoMerge = true;

		if (DoMerge)
			Prev.second = Blocks[i].second;
		else
			Merged.push_back(Blocks[i]);
		}
	Blocks = Merged;
	}

static void AppendCoreColumns(const MSA &Aln, uint Lo, uint Hi,
  vector<string> &NewRows, uint &CleanCol,
  vector<MSAPrepColMapEntry> &ColMap,
  const vector<vector<uint> > &OrigColToPos)
	{
	const uint SeqCount = Aln.GetSeqCount();
	for (uint Col = Lo; Col <= Hi; ++Col)
		{
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char c = Aln.GetChar(SeqIndex, Col);
			NewRows[SeqIndex] += c;
			if (!isgap(c))
				{
				MSAPrepColMapEntry E;
				E.CleanCol = CleanCol;
				E.OrigCol = Col;
				E.ResiduePos = OrigColToPos[SeqIndex][Col];
				E.SeqIndex = SeqIndex;
				ColMap.push_back(E);
				}
			}
		++CleanCol;
		}
	}

static uint MedianUint(vector<uint> Vals)
	{
	if (Vals.empty())
		return 0;
	sort(Vals.begin(), Vals.end());
	return Vals[SIZE(Vals)/2];
	}

static void AppendSpacerRegion(const MSA &Aln, uint Lo, uint Hi,
  const MSAPrepParams &Params,
  vector<string> &NewRows, uint &CleanCol,
  vector<MSAPrepColMapEntry> &ColMap,
  const vector<vector<uint> > &OrigColToPos,
  vector<uint> &OrigColsKept)
	{
	OrigColsKept.clear();
	if (Lo > Hi)
		return;

	const uint SeqCount = Aln.GetSeqCount();
	asserta(Params.SpacerMinCols >= 1);
	asserta(Params.SpacerMaxCols >= Params.SpacerMinCols);

	vector<uint> PresentLens;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint Len = 0;
		for (uint Col = Lo; Col <= Hi; ++Col)
			{
			char c = Aln.GetChar(SeqIndex, Col);
			if (!isgap(c))
				++Len;
			}
		if (Len > 0)
			PresentLens.push_back(Len);
		}

	uint MedianLen = MedianUint(PresentLens);
	uint W = Params.SpacerMinCols;
	if (MedianLen > 0)
		{
		W = MedianLen;
		if (W < Params.SpacerMinCols)
			W = Params.SpacerMinCols;
		if (W > Params.SpacerMaxCols)
			W = Params.SpacerMaxCols;
		}

	const uint N = Hi - Lo + 1;
	uint Start = Lo;
	uint Count = N;
	if (N > W)
		{
		Start = Lo + (N - W)/2;
		Count = W;
		}
	const uint End = Start + Count - 1;
	for (uint Col = Start; Col <= End; ++Col)
		OrigColsKept.push_back(Col);

	AppendCoreColumns(Aln, Start, End, NewRows, CleanCol,
	  ColMap, OrigColToPos);
	}

MSAPrepResult MSAPrep(const MSA &Input, const MSAPrepParams &Params)
	{
	MSAPrepResult Result;
	Result.OrigColCount = Input.GetColCount();
	const uint SeqCount = Input.GetSeqCount();
	if (SeqCount == 0 || Result.OrigColCount == 0)
		{
		Result.CleanedMSA.Copy(Input);
		Result.CleanColCount = Result.OrigColCount;
		return Result;
		}

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string UngappedSeq;
		Input.GetUngappedSeqStr(SeqIndex, UngappedSeq);
		Result.FullUngappedSeqs.push_back(UngappedSeq);
		}

	vector<vector<uint> > OrigColToPos(SeqCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		Input.GetColToPos(SeqIndex, OrigColToPos[SeqIndex]);

	bool UseLowerCase = false;
	if (Params.ColMode == MSAPrepColMode::LowerCaseInserts)
		UseLowerCase = true;
	else if (Params.ColMode == MSAPrepColMode::OccupancyOnly)
		UseLowerCase = false;
	else
		UseLowerCase = MSAHasMixedCaseCols(Input);

	vector<bool> IsCore;
	GetCoreColMask(Input, Params, UseLowerCase, IsCore);

	uint TrimLo = 0;
	uint TrimHi = 0;
	EndTrimRange(IsCore, TrimLo, TrimHi);

	vector<pair<uint, uint> > CoreBlocks;
	FindCoreBlocks(IsCore, TrimLo, TrimHi, Params.MinCoreBlockCols, CoreBlocks);
	MergeCoreBlocks(Input, CoreBlocks, Params);

	if (CoreBlocks.empty())
		{
		ProgressLog("msa_prep: no core blocks found, keeping end-trimmed occupancy columns\n");
		for (uint Col = TrimLo; Col <= TrimHi; ++Col)
			{
			if (IsCore[Col])
				CoreBlocks.push_back(make_pair(Col, Col));
			}
		if (CoreBlocks.empty())
			{
			Result.CleanedMSA.Copy(Input);
			Result.CleanColCount = Result.OrigColCount;
			return Result;
			}
		}

	vector<string> NewRows(SeqCount);
	uint CleanCol = 0;
	uint BlockId = 0;
	uint InsertId = 0;

	auto EmitTerminalDiscard = [&](uint Lo, uint Hi)
		{
		if (Lo > Hi)
			return;
		MSAPrepInsertRegion Ins;
		Ins.RegionId = InsertId++;
		Ins.OrigColLo = Lo;
		Ins.OrigColHi = Hi;
		Ins.m_Policy = MSAPrepInsertRegion::Policy::Discard;
		Ins.CleanWidth = 0;
		Ins.CleanColLo = UINT_MAX;
		Result.InsertRegions.push_back(Ins);
		};

	auto EmitInterBlockSpacer = [&](uint Lo, uint Hi)
		{
		if (Lo > Hi)
			return;
		MSAPrepInsertRegion Ins;
		Ins.RegionId = InsertId++;
		Ins.OrigColLo = Lo;
		Ins.OrigColHi = Hi;
		Ins.m_Policy = MSAPrepInsertRegion::Policy::Spacer;
		Ins.CleanColLo = CleanCol;
		AppendSpacerRegion(Input, Lo, Hi, Params,
		  NewRows, CleanCol, Result.ColMap, OrigColToPos, Ins.OrigCols);
		Ins.CleanWidth = CleanCol - Ins.CleanColLo;
		asserta(SIZE(Ins.OrigCols) == Ins.CleanWidth);
		Result.InsertRegions.push_back(Ins);
		};

	if (CoreBlocks[0].first > 0)
		EmitTerminalDiscard(0, CoreBlocks[0].first - 1);

	for (uint bi = 0; bi < SIZE(CoreBlocks); ++bi)
		{
		uint BLo = CoreBlocks[bi].first;
		uint BHi = CoreBlocks[bi].second;

		MSAPrepBlock Block;
		Block.BlockId = BlockId++;
		Block.OrigColLo = BLo;
		Block.OrigColHi = BHi;
		Block.CleanColLo = CleanCol;
		AppendCoreColumns(Input, BLo, BHi, NewRows, CleanCol,
		  Result.ColMap, OrigColToPos);
		Block.CleanColHi = CleanCol - 1;
		Result.Blocks.push_back(Block);

		if (bi + 1 < SIZE(CoreBlocks))
			{
			uint GapLo = BHi + 1;
			uint GapHi = CoreBlocks[bi + 1].first - 1;
			EmitInterBlockSpacer(GapLo, GapHi);
			}
		}

	if (CoreBlocks.back().second + 1 < Result.OrigColCount)
		EmitTerminalDiscard(CoreBlocks.back().second + 1, Result.OrigColCount - 1);

	vector<string> Labels;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		Labels.push_back(string(Input.GetLabel(SeqIndex)));

	Result.CleanedMSA.FromStrings2(Labels, NewRows);
	Result.CleanColCount = CleanCol;
	for (uint i = 0; i < SIZE(Result.ColMap); ++i)
		asserta(Result.ColMap[i].OrigCol != UINT_MAX);

	ProgressLog("msa_prep: %u -> %u cols, %u blocks, %u insert regions\n",
	  Result.OrigColCount, Result.CleanColCount,
	  SIZE(Result.Blocks), SIZE(Result.InsertRegions));

	return Result;
	}

void WriteMSAPrepBlocksFile(const string &FileName,
  const MSAPrepResult &Result, const MSAPrepParams &Params)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "msa_prep_blocks\t1\n");
	fprintf(f, "params\tmax_gap_fract\t%.4g\tmin_core_block_cols\t%u\t"
	  "spacer_min_cols\t%u\tspacer_max_cols\t%u\t"
	  "spacer_min_occ\t%.4g\tspacer_max_occ\t%.4g\t"
	  "merge_gap_cols\t%u\n",
	  Params.MaxGapFract, Params.MinCoreBlockCols,
	  Params.SpacerMinCols, Params.SpacerMaxCols,
	  Params.SpacerMinOcc, Params.SpacerMaxOcc,
	  Params.MergeGapCols);
	fprintf(f, "orig_cols\t%u\n", Result.OrigColCount);
	fprintf(f, "clean_cols\t%u\n", Result.CleanColCount);
	for (uint i = 0; i < SIZE(Result.Blocks); ++i)
		{
		const MSAPrepBlock &B = Result.Blocks[i];
		fprintf(f, "block\t%u\torig_lo\t%u\torig_hi\t%u\tclean_lo\t%u\tclean_hi\t%u\n",
		  B.BlockId, B.OrigColLo, B.OrigColHi, B.CleanColLo, B.CleanColHi);
		}
	for (uint i = 0; i < SIZE(Result.InsertRegions); ++i)
		{
		const MSAPrepInsertRegion &Ins = Result.InsertRegions[i];
		const char *PolicyStr = "discard";
		if (Ins.m_Policy == MSAPrepInsertRegion::Policy::Compress)
			PolicyStr = "compress";
		else if (Ins.m_Policy == MSAPrepInsertRegion::Policy::Spacer)
			PolicyStr = "spacer";
		if (Ins.CleanColLo == UINT_MAX)
			fprintf(f, "insert\t%u\torig_lo\t%u\torig_hi\t%u\tpolicy\t%s\t"
			  "clean_cols\t%u\tclean_lo\t-\n",
			  Ins.RegionId, Ins.OrigColLo, Ins.OrigColHi, PolicyStr,
			  Ins.CleanWidth);
		else
			fprintf(f, "insert\t%u\torig_lo\t%u\torig_hi\t%u\tpolicy\t%s\t"
			  "clean_cols\t%u\tclean_lo\t%u\n",
			  Ins.RegionId, Ins.OrigColLo, Ins.OrigColHi, PolicyStr,
			  Ins.CleanWidth, Ins.CleanColLo);
		}
	CloseStdioFile(f);
	}

void WriteMSAPrepMapFile(const string &FileName,
  const MSAPrepResult &Result)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "msa_prep_map\t1\n");
	fprintf(f, "orig_cols\t%u\n", Result.OrigColCount);
	fprintf(f, "seed_cols\t%u\n", Result.CleanColCount);

	uint BlockIdx = 0;
	uint InsIdx = 0;
	uint ListedSeedCols = 0;
	uint SpacerId = 0;
	const uint BlockCount = SIZE(Result.Blocks);
	const uint InsCount = SIZE(Result.InsertRegions);
	while (BlockIdx < BlockCount || InsIdx < InsCount)
		{
		const bool HaveBlock = (BlockIdx < BlockCount);
		const bool HaveIns = (InsIdx < InsCount);
		const bool EmitBlock = HaveBlock &&
		  (!HaveIns || Result.Blocks[BlockIdx].OrigColLo <=
		    Result.InsertRegions[InsIdx].OrigColLo);
		if (EmitBlock)
			{
			const MSAPrepBlock &B = Result.Blocks[BlockIdx++];
			fprintf(f, "block\t%u\torig_lo\t%u\torig_hi\t%u\tseed_lo\t%u\tseed_hi\t%u\n",
			  B.BlockId, B.OrigColLo, B.OrigColHi, B.CleanColLo, B.CleanColHi);
			asserta(B.OrigColHi >= B.OrigColLo);
			ListedSeedCols += B.OrigColHi - B.OrigColLo + 1;
			}
		else
			{
			const MSAPrepInsertRegion &Ins = Result.InsertRegions[InsIdx++];
			if (Ins.m_Policy == MSAPrepInsertRegion::Policy::Discard ||
			  Ins.CleanWidth == 0 || Ins.CleanColLo == UINT_MAX)
				{
				fprintf(f, "discard\torig_lo\t%u\torig_hi\t%u\n",
				  Ins.OrigColLo, Ins.OrigColHi);
				}
			else
				{
				uint SeedHi = Ins.CleanColLo + Ins.CleanWidth - 1;
				fprintf(f, "spacer\t%u\torig_lo\t%u\torig_hi\t%u\t"
				  "seed_lo\t%u\tseed_hi\t%u\torig_cols\t",
				  SpacerId++, Ins.OrigColLo, Ins.OrigColHi,
				  Ins.CleanColLo, SeedHi);
				for (uint i = 0; i < SIZE(Ins.OrigCols); ++i)
					{
					if (i > 0)
						fprintf(f, ",");
					fprintf(f, "%u", Ins.OrigCols[i]);
					}
				fprintf(f, "\n");
				ListedSeedCols += SIZE(Ins.OrigCols);
				}
			}
		}
	asserta(ListedSeedCols == Result.CleanColCount);
	CloseStdioFile(f);
	}

MSAPrepParams GetMSAPrepParamsFromOpts()
	{
	MSAPrepParams Params;
	if (optset_occupancy_only)
		Params.ColMode = MSAPrepColMode::OccupancyOnly;
	else if (optset_lower_case_inserts)
		Params.ColMode = MSAPrepColMode::LowerCaseInserts;
	else
		Params.ColMode = MSAPrepColMode::Auto;

	if (optset_max_gap_fract)
		Params.MaxGapFract = opt(max_gap_fract);
	if (optset_min_core_block_cols)
		Params.MinCoreBlockCols = opt(min_core_block_cols);
	if (optset_insert_minority_frac)
		Params.InsertMinorityFrac = opt(insert_minority_frac);
	if (optset_insert_discard_len)
		Params.InsertDiscardLen = opt(insert_discard_len);
	if (optset_insert_cap_len)
		Params.InsertCapLen = opt(insert_cap_len);
	if (optset_merge_gap_cols)
		Params.MergeGapCols = opt(merge_gap_cols);
	if (optset_spacer_min_cols)
		Params.SpacerMinCols = opt(spacer_min_cols);
	if (optset_spacer_max_cols)
		Params.SpacerMaxCols = opt(spacer_max_cols);
	if (optset_spacer_min_occ)
		Params.SpacerMinOcc = opt(spacer_min_occ);
	if (optset_spacer_max_occ)
		Params.SpacerMaxOcc = opt(spacer_max_occ);
	if (Params.SpacerMinCols > Params.SpacerMaxCols)
		Die("spacer_min_cols > spacer_max_cols");
	if (Params.SpacerMinOcc > Params.SpacerMaxOcc)
		Die("spacer_min_occ > spacer_max_occ");
	return Params;
	}

static const char *g_MSAPrepBlockColors_JalView[] =
	{
	"0040A0",
	"A04000",
	"008040",
	"804000",
	"400080",
	"408000",
	"004080",
	"800040",
	"806000",
	"006080",
	"608000",
	"800060"
	};
static const uint g_MSAPrepBlockColorCount = 12;

void WriteMSAPrepBlocksJalView(const string &FileName,
  const MSA &OrigAln, const MSAPrepResult &Result)
	{
	if (FileName.empty())
		return;

	const uint SeqCount = OrigAln.GetSeqCount();
	FILE *fOut = CreateStdioFile(FileName);
	for (uint i = 0; i < g_MSAPrepBlockColorCount; ++i)
		fprintf(fOut, "Block%u\t%s\n", i, g_MSAPrepBlockColors_JalView[i]);

	fprintf(fOut, "STARTGROUP\tMSAPrep_CoreBlocks\n");
	for (uint bi = 0; bi < SIZE(Result.Blocks); ++bi)
		{
		const MSAPrepBlock &B = Result.Blocks[bi];
		uint ColorIdx = bi % g_MSAPrepBlockColorCount;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			vector<uint> ColToPos;
			OrigAln.GetColToPos(SeqIndex, ColToPos);
			uint SeqLo = UINT_MAX;
			uint SeqHi = UINT_MAX;
			for (uint Col = B.OrigColLo; Col <= B.OrigColHi; ++Col)
				{
				if (OrigAln.IsGap(SeqIndex, Col))
					continue;
				uint Pos = ColToPos[Col];
				if (SeqLo == UINT_MAX)
					SeqLo = Pos;
				SeqHi = Pos;
				}
			if (SeqLo == UINT_MAX)
				continue;
			string Label;
			OrigAln.GetSeqLabel(SeqIndex, Label);
			fprintf(fOut, "-\t%s\t%u\t%u\t%u\tBlock%u\n",
			  Label.c_str(), SeqIndex, SeqLo + 1, SeqHi + 1, ColorIdx);
			}
		}
	fprintf(fOut, "ENDGROUP\tMSAPrep_CoreBlocks\n");
	CloseStdioFile(fOut);
	}

void MSAToMultiSequence(const MSA &Aln, MultiSequence &MS)
	{
	vector<string> Labels;
	vector<string> Rows;
	const uint SeqCount = Aln.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		Labels.push_back(string(Aln.GetLabel(SeqIndex)));
		string Row;
		Aln.GetRowStr(SeqIndex, Row);
		Rows.push_back(Row);
		}
	MS.FromStrings(Labels, Rows);
	}
