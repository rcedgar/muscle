#include "muscle.h"
#include "masm.h"
#include "xdpmem.h"

static const uint BIN_COUNT = 100;
static const uint EMPTY_QUERY_TRIES = 100;
static const uint BLOCK_PERM_TRIES = 1000;

float SWFast_MASM(XDPMem &Mem, const MASM &A, const vector<vector<byte> > &B,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);
float NWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, uint &Loj, string &Path);

struct CalibBlock
	{
	uint SeedLo = 0;
	uint SeedHi = 0;
	};

static float RandUnit()
	{
	return float(randu32())*(1.0f/4294967296.0f);
	}

static float SumFreqs(const vector<float> &Freqs)
	{
	float Sum = 0;
	for (uint i = 0; i < SIZE(Freqs); ++i)
		Sum += Freqs[i];
	return Sum;
	}

static byte SampleLetter(const vector<float> &Freqs, float Occ)
	{
	asserta(Occ > 0);
	const uint AlphaSize = SIZE(Freqs);
	float v = RandUnit()*Occ;
	float Cum = 0;
	byte Last = 0;
	bool Any = false;
	for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		{
		float f = Freqs[Letter];
		if (f <= 0)
			continue;
		Last = (byte) Letter;
		Any = true;
		Cum += f;
		if (v < Cum)
			return (byte) Letter;
		}
	if (!Any)
		Die("calibrate_masm: empty feature freqs");
	return Last;
	}

static bool SampleMegaPos(const MASM &M, uint Col, vector<byte> &Pos)
	{
	asserta(Col < M.m_ColCount);
	asserta(M.m_AAFeatureIdx < M.m_FeatureCount);
	const MASMCol &MC = M.GetCol(Col);
	const vector<float> &AAFreqs = MC.m_FreqsVec[M.m_AAFeatureIdx];
	const float Occ = SumFreqs(AAFreqs);
	if (RandUnit() >= Occ)
		return false;

	const uint FeatureCount = M.m_FeatureCount;
	Pos.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		const vector<float> &Freqs = MC.m_FreqsVec[FeatureIdx];
		float FeatOcc = SumFreqs(Freqs);
		if (FeatOcc <= 0)
			Die("calibrate_masm: feature %u col %u has zero occupancy",
			  FeatureIdx, Col);
		Pos[FeatureIdx] = SampleLetter(Freqs, FeatOcc);
		}
	return true;
	}

static void SampleQueryFromCols(const MASM &M, const vector<uint> &ColOrder,
  vector<vector<byte> > &Prof)
	{
	Prof.clear();
	for (uint i = 0; i < SIZE(ColOrder); ++i)
		{
		vector<byte> Pos;
		if (SampleMegaPos(M, ColOrder[i], Pos))
			Prof.push_back(Pos);
		}
	}

static void SequentialCols(uint ColCount, vector<uint> &Cols)
	{
	Cols.resize(ColCount);
	for (uint i = 0; i < ColCount; ++i)
		Cols[i] = i;
	}

static bool TrySampleQuery(const MASM &M, const vector<uint> &ColOrder,
  vector<vector<byte> > &Prof)
	{
	for (uint Try = 0; Try < EMPTY_QUERY_TRIES; ++Try)
		{
		SampleQueryFromCols(M, ColOrder, Prof);
		if (!Prof.empty())
			return true;
		}
	return false;
	}

static float ScoreQuery(XDPMem &Mem, const MASM &Target,
  const vector<vector<byte> > &Prof, bool Local)
	{
	string Path;
	uint Loj;
	if (Local)
		{
		uint Loi, Leni, Lenj;
		return SWFast_MASM(Mem, Target, Prof, Loi, Loj, Leni, Lenj, Path);
		}
	return NWFast_MASM_MegaProf(Mem, Target, Prof, Loj, Path);
	}

static void AssertSameFeatures(const MASM &A, const MASM &B)
	{
	if (A.m_FeatureCount != B.m_FeatureCount)
		Die("calibrate_masm -decoy: feature count %u != target %u",
		  B.m_FeatureCount, A.m_FeatureCount);
	for (uint i = 0; i < A.m_FeatureCount; ++i)
		{
		if (A.m_FeatureNames[i] != B.m_FeatureNames[i])
			Die("calibrate_masm -decoy: feature %u name '%s' != target '%s'",
			  i, B.m_FeatureNames[i].c_str(), A.m_FeatureNames[i].c_str());
		if (A.m_AlphaSizes[i] != B.m_AlphaSizes[i])
			Die("calibrate_masm -decoy: feature %u alpha %u != target %u",
			  i, B.m_AlphaSizes[i], A.m_AlphaSizes[i]);
		}
	}

static void ParseMapFile(const string &FileName, uint ColCount,
  vector<CalibBlock> &Blocks, vector<uint> &InsertPool,
  vector<uint> &SpacerWidths)
	{
	Blocks.clear();
	InsertPool.clear();
	SpacerWidths.clear();

	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);
	if (Lines.empty())
		Die("calibrate_masm: empty map file '%s'", FileName.c_str());

	vector<string> Fields;
	for (uint LineIdx = 0; LineIdx < SIZE(Lines); ++LineIdx)
		{
		const string &Line = Lines[LineIdx];
		if (Line.empty())
			continue;
		Split(Line, Fields, '\t');
		if (Fields.empty())
			continue;
		const string &Key = Fields[0];
		if (Key == "msa_prep_map" || Key == "orig_cols" || Key == "seed_cols" ||
		  Key == "discard")
			continue;
		if (Key == "block")
			{
			if (SIZE(Fields) < 10)
				Die("calibrate_masm: bad block line in '%s'", FileName.c_str());
			asserta(Fields[2] == "orig_lo");
			asserta(Fields[4] == "orig_hi");
			asserta(Fields[6] == "seed_lo");
			asserta(Fields[8] == "seed_hi");
			CalibBlock B;
			B.SeedLo = StrToUint(Fields[7]);
			B.SeedHi = StrToUint(Fields[9]);
			if (B.SeedHi < B.SeedLo)
				Die("calibrate_masm: block seed_hi < seed_lo in '%s'",
				  FileName.c_str());
			if (B.SeedHi >= ColCount)
				Die("calibrate_masm: block seed_hi %u >= MASM cols %u",
				  B.SeedHi, ColCount);
			Blocks.push_back(B);
			continue;
			}
		if (Key == "spacer")
			{
			if (SIZE(Fields) < 10)
				Die("calibrate_masm: bad spacer line in '%s'", FileName.c_str());
			asserta(Fields[6] == "seed_lo");
			asserta(Fields[8] == "seed_hi");
			uint SeedLo = StrToUint(Fields[7]);
			uint SeedHi = StrToUint(Fields[9]);
			if (SeedHi < SeedLo)
				Die("calibrate_masm: spacer seed_hi < seed_lo in '%s'",
				  FileName.c_str());
			if (SeedHi >= ColCount)
				Die("calibrate_masm: spacer seed_hi %u >= MASM cols %u",
				  SeedHi, ColCount);
			SpacerWidths.push_back(SeedHi - SeedLo + 1);
			for (uint Col = SeedLo; Col <= SeedHi; ++Col)
				InsertPool.push_back(Col);
			continue;
			}
		}

	if (Blocks.empty())
		Die("calibrate_masm -denovo: no block rows in '%s'", FileName.c_str());

	const uint BlockCount = SIZE(Blocks);
	for (uint i = 0; i < BlockCount; ++i)
		{
		uint Best = i;
		for (uint j = i + 1; j < BlockCount; ++j)
			if (Blocks[j].SeedLo < Blocks[Best].SeedLo)
				Best = j;
		if (Best != i)
			{
			CalibBlock Tmp = Blocks[i];
			Blocks[i] = Blocks[Best];
			Blocks[Best] = Tmp;
			}
		}
	}

static bool PermHasOrigAdjacent(const vector<uint> &Perm)
	{
	const uint K = SIZE(Perm);
	for (uint i = 0; i + 1 < K; ++i)
		if (Perm[i] + 1 == Perm[i + 1])
			return true;
	return false;
	}

static void PermuteBlocks(uint K, vector<uint> &Perm)
	{
	Perm.resize(K);
	for (uint i = 0; i < K; ++i)
		Perm[i] = i;
	if (K < 2)
		return;

	for (uint Try = 0; Try < BLOCK_PERM_TRIES; ++Try)
		{
		Shuffle(Perm);
		if (!PermHasOrigAdjacent(Perm))
			return;
		}
	Die("calibrate_masm -denovo: failed to permute %u core blocks", K);
	}

static void MakeDenovoColOrder(const vector<CalibBlock> &Blocks,
  const vector<uint> &InsertPool, const vector<uint> &SpacerWidths,
  vector<uint> &ColOrder)
	{
	const uint K = SIZE(Blocks);
	vector<uint> Perm;
	PermuteBlocks(K, Perm);

	ColOrder.clear();
	for (uint pi = 0; pi < K; ++pi)
		{
		if (pi > 0 && !SpacerWidths.empty() && !InsertPool.empty())
			{
			uint W = SpacerWidths[randu32()%SIZE(SpacerWidths)];
			for (uint w = 0; w < W; ++w)
				ColOrder.push_back(InsertPool[randu32()%SIZE(InsertPool)]);
			}
		const CalibBlock &B = Blocks[Perm[pi]];
		for (uint Col = B.SeedLo; Col <= B.SeedHi; ++Col)
			ColOrder.push_back(Col);
		}
	}

static void MakeShatterColOrder(uint ColCount, uint Min, uint Max,
  vector<uint> &ColOrder)
	{
	if (ColCount < 2)
		Die("calibrate_masm -shatter: need at least 2 columns");
	asserta(Min >= 1 && Min <= Max);

	vector<uint> SegLo;
	vector<uint> SegLen;
	uint Pos = 0;
	while (Pos < ColCount)
		{
		uint Rem = ColCount - Pos;
		uint Span = Max - Min + 1;
		uint L = Min + randu32()%Span;
		if (L > Rem)
			L = Rem;
		SegLo.push_back(Pos);
		SegLen.push_back(L);
		Pos += L;
		}

	const uint SegCount = SIZE(SegLo);
	vector<uint> Order(SegCount);
	for (uint i = 0; i < SegCount; ++i)
		Order[i] = i;
	Shuffle(Order);

	ColOrder.clear();
	for (uint k = 0; k < SegCount; ++k)
		{
		uint i = Order[k];
		for (uint d = 0; d < SegLen[i]; ++d)
			ColOrder.push_back(SegLo[i] + d);
		}
	asserta(SIZE(ColOrder) == ColCount);
	}

static void ScoreStats(const vector<float> &Scores,
  float &Min, float &Max, float &Mean)
	{
	const uint N = SIZE(Scores);
	asserta(N > 0);
	Min = Scores[0];
	Max = Scores[0];
	double Sum = 0;
	for (uint i = 0; i < N; ++i)
		{
		float x = Scores[i];
		if (x < Min)
			Min = x;
		if (x > Max)
			Max = x;
		Sum += x;
		}
	Mean = float(Sum/N);
	}

static double RoundToPowerOf10(double P)
	{
	if (P <= 0)
		Die("calibrate: P-value at lowtp is not positive");
	if (P >= 1)
		return 1;
	double e = log10(P);
	return pow(10.0, round(e));
	}

static void WriteHistTSV(const string &FileName,
  const vector<float> &TPScores, const vector<float> &FPScores,
  double &Slope, double &Intercept, double &Hi,
  double &LowTP, double &Cutoff)
	{
	const uint NTP = SIZE(TPScores);
	const uint NFP = SIZE(FPScores);
	asserta(NTP > 0 && NFP > 0);

	float HistLo = TPScores[0];
	float HistHi = TPScores[0];
	for (uint i = 0; i < NTP; ++i)
		{
		if (TPScores[i] < HistLo)
			HistLo = TPScores[i];
		if (TPScores[i] > HistHi)
			HistHi = TPScores[i];
		}
	for (uint i = 0; i < NFP; ++i)
		{
		if (FPScores[i] < HistLo)
			HistLo = FPScores[i];
		if (FPScores[i] > HistHi)
			HistHi = FPScores[i];
		}
	if (HistLo == HistHi)
		{
		HistLo -= 1e-3f;
		HistHi += 1e-3f;
		}

	const float Width = (HistHi - HistLo)/float(BIN_COUNT);
	vector<uint> TPCounts(BIN_COUNT, 0);
	vector<uint> FPCounts(BIN_COUNT, 0);

	for (uint i = 0; i < NTP; ++i)
		{
		int Bin = int((TPScores[i] - HistLo)/Width);
		if (Bin < 0)
			Bin = 0;
		if (Bin >= int(BIN_COUNT))
			Bin = int(BIN_COUNT) - 1;
		TPCounts[(uint) Bin] += 1;
		}
	for (uint i = 0; i < NFP; ++i)
		{
		int Bin = int((FPScores[i] - HistLo)/Width);
		if (Bin < 0)
			Bin = 0;
		if (Bin >= int(BIN_COUNT))
			Bin = int(BIN_COUNT) - 1;
		FPCounts[(uint) Bin] += 1;
		}

	uint TailHi = UINT_MAX;
	uint TailLo = UINT_MAX;
	for (uint Bin = 0; Bin < BIN_COUNT; ++Bin)
		{
		if (FPCounts[Bin] >= 100)
			TailHi = Bin;
		if (FPCounts[Bin] >= 10)
			TailLo = Bin;
		}
	if (TailHi == UINT_MAX)
		Die("calibrate_masm: no FP bin with >= 100 scores");
	if (TailLo == UINT_MAX || TailLo <= TailHi)
		Die("calibrate_masm: FP tail bins invalid (hi %u, lo %u)",
		  TailHi, TailLo);

	double SumX = 0;
	double SumY = 0;
	double SumXX = 0;
	double SumXY = 0;
	uint NFit = 0;
	for (uint Bin = TailHi; Bin <= TailLo; ++Bin)
		{
		uint y = FPCounts[Bin];
		if (y == 0)
			continue;
		double x = double(HistLo) + (double(Bin) + 0.5)*double(Width);
		double ly = log(double(y));
		SumX += x;
		SumY += ly;
		SumXX += x*x;
		SumXY += x*ly;
		++NFit;
		}
	if (NFit < 2)
		Die("calibrate_masm: need >= 2 positive FP bins in tail [%u, %u]",
		  TailHi, TailLo);
	double Det = double(NFit)*SumXX - SumX*SumX;
	if (Det == 0)
		Die("calibrate_masm: log-linear fit degenerate");
	double b = (double(NFit)*SumXY - SumX*SumY)/Det;
	double a = (SumY - b*SumX)/double(NFit);
	Slope = b;
	Intercept = a;
	Hi = double(HistLo) + (double(TailHi) + 0.5)*double(Width);

	uint LowTPBin = UINT_MAX;
	for (uint Bin = 0; Bin < BIN_COUNT; ++Bin)
		{
		if (TPCounts[Bin] >= 10)
			{
			LowTPBin = Bin;
			break;
			}
		}
	if (LowTPBin == UINT_MAX)
		Die("calibrate_masm: no TP bin with >= 10 scores");
	LowTP = double(HistLo) + (double(LowTPBin) + 0.5)*double(Width);
	double P = exp(a + b*LowTP)/double(NFP);
	if (P > 1)
		P = 1;
	Cutoff = RoundToPowerOf10(P);
	ProgressLog("FP loglin  bins %u..%u  n %u  a %.4g  b %.4g  hi %.4g\n",
	  TailHi, TailLo, NFit, a, b, Hi);
	ProgressLog("TP lowtp  bin %u  score %.4g  P %.4g  cutoff %.4g\n",
	  LowTPBin, LowTP, P, Cutoff);

	if (!FileName.empty())
		{
		FILE *f = CreateStdioFile(FileName);
		fprintf(f, "bin\tlo\thi\ttp\tfp\tloglin\n");
		for (uint Bin = 0; Bin < BIN_COUNT; ++Bin)
			{
			float BinLo = HistLo + float(Bin)*Width;
			float BinHi = HistLo + float(Bin + 1)*Width;
			fprintf(f, "%u\t%.6g\t%.6g\t%u\t%u\t",
			  Bin, BinLo, BinHi, TPCounts[Bin], FPCounts[Bin]);
			if (Bin < TailHi)
				fprintf(f, "-\n");
			else
				{
				double x = double(HistLo) + (double(Bin) + 0.5)*double(Width);
				double Fit = exp(a + b*x);
				if (Fit < 0)
					Fit = 0;
				fprintf(f, "%.4g\n", Fit);
				}
			}
		CloseStdioFile(f);
		}
	}

void CalibrateMASM(MASM &Target, bool Local, uint N,
  bool Denovo, const string &MapFN,
  bool Decoy, const string &DecoyFN,
  bool Shatter, uint ShatterMin, uint ShatterMax,
  const string &HistTSVFN)
	{
	uint FPModeCount = 0;
	if (Denovo)
		++FPModeCount;
	if (Decoy)
		++FPModeCount;
	if (Shatter)
		++FPModeCount;
	if (FPModeCount != 1)
		Die("calibrate: require exactly one of denovo, decoy or shatter");
	if (N == 0)
		Die("calibrate: -n must be > 0");
	if (Target.m_ColCount == 0)
		Die("calibrate: MASM has 0 columns");
	if (Target.m_AAFeatureIdx == UINT_MAX)
		Die("calibrate: MASM has no AA feature");
	if (Target.m_FeatureCount == 0)
		Die("calibrate: MASM has no features");
	if (Shatter)
		{
		if (ShatterMin < 1 || ShatterMin > ShatterMax)
			Die("calibrate -shatter: bad range %u..%u",
			  ShatterMin, ShatterMax);
		if (Target.m_ColCount < 2)
			Die("calibrate -shatter: need at least 2 columns");
		}

	MASM DecoyMASM;
	const MASM *SampleMASM = &Target;
	vector<CalibBlock> Blocks;
	vector<uint> InsertPool;
	vector<uint> SpacerWidths;
	if (Denovo)
		{
		if (MapFN.empty() || !StdioFileExists(MapFN))
			Die("calibrate -denovo: map file not found '%s'",
			  MapFN.c_str());
		ParseMapFile(MapFN, Target.m_ColCount, Blocks, InsertPool,
		  SpacerWidths);
		}
	else if (Decoy)
		{
		if (DecoyFN.empty())
			Die("calibrate: -decoy requires a MASM file");
		DecoyMASM.FromFile(DecoyFN);
		AssertSameFeatures(Target, DecoyMASM);
		if (DecoyMASM.m_ColCount == 0)
			Die("calibrate: decoy MASM has 0 columns");
		if (DecoyMASM.m_AAFeatureIdx == UINT_MAX)
			Die("calibrate: decoy MASM has no AA feature");
		SampleMASM = &DecoyMASM;
		}

	vector<uint> TPCols;
	SequentialCols(Target.m_ColCount, TPCols);
	vector<uint> DecoyCols;
	if (Decoy)
		SequentialCols(SampleMASM->m_ColCount, DecoyCols);

	XDPMem Mem;
	vector<float> TPScores;
	vector<float> FPScores;
	TPScores.reserve(N);
	FPScores.reserve(N);

	vector<vector<byte> > Prof;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Calibrating TP");
		if (!TrySampleQuery(Target, TPCols, Prof))
			Die("calibrate: empty TP query after %u tries",
			  EMPTY_QUERY_TRIES);
		TPScores.push_back(ScoreQuery(Mem, Target, Prof, Local));
		}

	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Calibrating FP");
		vector<uint> FPCols;
		if (Denovo)
			MakeDenovoColOrder(Blocks, InsertPool, SpacerWidths, FPCols);
		else if (Shatter)
			MakeShatterColOrder(Target.m_ColCount, ShatterMin, ShatterMax,
			  FPCols);
		else
			FPCols = DecoyCols;
		if (!TrySampleQuery(*SampleMASM, FPCols, Prof))
			Die("calibrate: empty FP query after %u tries",
			  EMPTY_QUERY_TRIES);
		FPScores.push_back(ScoreQuery(Mem, Target, Prof, Local));
		}

	float TPMin, TPMax, TPMean;
	float FPMin, FPMax, FPMean;
	ScoreStats(TPScores, TPMin, TPMax, TPMean);
	ScoreStats(FPScores, FPMin, FPMax, FPMean);
	ProgressLog("TP  n %u  min %.4g  mean %.4g  max %.4g\n",
	  N, TPMin, TPMean, TPMax);
	ProgressLog("FP  n %u  min %.4g  mean %.4g  max %.4g\n",
	  N, FPMin, FPMean, FPMax);

	double Slope = 0;
	double Intercept = 0;
	double CalibHi = 0;
	double LowTP = 0;
	double Cutoff = 1e-3;
	WriteHistTSV(HistTSVFN, TPScores, FPScores, Slope, Intercept, CalibHi,
	  LowTP, Cutoff);
	Target.m_HasCalibrate = true;
	Target.m_CalibSlope = Slope;
	Target.m_CalibIntercept = Intercept;
	Target.m_CalibHi = CalibHi;
	Target.m_CalibSamples = N;
	Target.m_CalibLowTP = LowTP;
	Target.m_CalibCutoff = Cutoff;
	if (!HistTSVFN.empty())
		ProgressLog("Wrote %s\n", HistTSVFN.c_str());
	}

void cmd_calibrate_masm()
	{
	if (opt(local) == opt(global))
		Die("calibrate_masm: require exactly one of -local or -global");
	uint FPModeCount = 0;
	if (opt(denovo))
		++FPModeCount;
	if (optset_decoy)
		++FPModeCount;
	if (opt(shatter))
		++FPModeCount;
	if (FPModeCount != 1)
		Die("calibrate_masm: require exactly one of -denovo, -decoy or -shatter");
	if (optset_decoy && opt(decoy).empty())
		Die("calibrate_masm: -decoy requires a MASM file");
	if (!optset_output)
		Die("calibrate_masm: -output required");
	if (!optset_tsvout)
		Die("calibrate_masm: -tsvout required");
	if (g_Arg1.empty())
		Die("calibrate_masm: missing target MASM");

	const bool Local = opt(local);
	const uint N = optset_n ? opt(n) : 10000;
	if (N == 0)
		Die("calibrate_masm: -n must be > 0");

	uint ShatterMin = 5;
	uint ShatterMax = 15;
	if (opt(shatter))
		{
		if (optset_shatter_min)
			ShatterMin = opt(shatter_min);
		if (optset_shatter_max)
			ShatterMax = opt(shatter_max);
		}

	string MapFN;
	if (opt(denovo))
		{
		if (optset_map)
			MapFN = opt(map);
		else
			MapFN = g_Arg1 + ".map";
		}

	MASM Target;
	Target.FromFile(g_Arg1);
	CalibrateMASM(Target, Local, N,
	  opt(denovo), MapFN,
	  optset_decoy, optset_decoy ? opt(decoy) : "",
	  opt(shatter), ShatterMin, ShatterMax,
	  opt(tsvout));
	Target.ToFile(opt(output));
	ProgressLog("Wrote %s\n", opt(output).c_str());
	}

void cmd_strumm_calibrate()
	{
	cmd_calibrate_masm();
	}
