#include "muscle.h"
#include "masm.h"
#include "xdpmem.h"
#include "omplock.h"

float SWFast_MASM(XDPMem &Mem, const MASM &A, const vector<vector<byte> > &B,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);
float NWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, uint &Loj, string &Path);
XDPMem &GetDPMem();

static void AssertMASMMatchesMega(const MASM &M)
	{
	if (M.m_FeatureCount != Mega::GetFeatureCount())
		Die("masm_search: MASM feature count %u != STRUCTS %u",
		  M.m_FeatureCount, Mega::GetFeatureCount());
	for (uint i = 0; i < M.m_FeatureCount; ++i)
		{
		if (M.m_FeatureNames[i] != Mega::GetFeatureName(i))
			Die("masm_search: MASM feature %u name '%s' != STRUCTS '%s'",
			  i, M.m_FeatureNames[i].c_str(), Mega::GetFeatureName(i).c_str());
		if (M.m_AlphaSizes[i] != Mega::GetAlphaSize(i))
			Die("masm_search: MASM feature %u alpha %u != STRUCTS %u",
			  i, M.m_AlphaSizes[i], Mega::GetAlphaSize(i));
		}
	}

static char A3MUpper(char c)
	{
	if (isalpha(c))
		return toupper(c);
	return c;
	}

static char A3MLower(char c)
	{
	if (isalpha(c))
		return tolower(c);
	return c;
	}

static void PathToA3M(const string &Seq, uint LA, uint Loi, uint Loj,
  const string &Path, string &A3M)
	{
	const uint LQ = SIZE(Seq);
	A3M.clear();
	if (Path.empty())
		{
		for (uint i = 0; i < LQ; ++i)
			A3M += A3MLower(Seq[i]);
		A3M.append(LA, '-');
		return;
		}

	asserta(Loj <= LQ);
	asserta(Loi <= LA);
	for (uint i = 0; i < Loj; ++i)
		A3M += A3MLower(Seq[i]);
	for (uint i = 0; i < Loi; ++i)
		A3M += '-';

	uint PosA = Loi;
	uint PosB = Loj;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			asserta(PosA < LA);
			asserta(PosB < LQ);
			A3M += A3MUpper(Seq[PosB]);
			++PosA;
			++PosB;
			}
		else if (c == 'D')
			{
			asserta(PosA < LA);
			A3M += '-';
			++PosA;
			}
		else if (c == 'I')
			{
			asserta(PosB < LQ);
			A3M += A3MLower(Seq[PosB]);
			++PosB;
			}
		else
			Die("masm_search: bad path char '%c'", c);
		}

	while (PosA < LA)
		{
		A3M += '-';
		++PosA;
		}
	while (PosB < LQ)
		{
		A3M += A3MLower(Seq[PosB]);
		++PosB;
		}
	}

void cmd_masm_search()
	{
	if (g_Arg1.empty())
		Die("masm_search: missing STRUCTS");
	if (!optset_masm)
		Die("masm_search: -masm required");
	if (opt(local) && opt(global))
		Die("masm_search: specify at most one of -local or -global");
	const bool Local = opt(local);
	if (optset_global)
		(void) opt(global);

	const string &StructsFN = g_Arg1;
	Mega::RejectLegacyMega(StructsFN);
	Mega::FromStructs(StructsFN);

	MASM M;
	M.FromFile(opt(masm));
	if (M.m_ColCount == 0)
		Die("masm_search: MASM has 0 columns");
	AssertMASMMatchesMega(M);

	const uint QueryCount = Mega::GetProfileCount();
	if (QueryCount == 0)
		Die("masm_search: no queries in STRUCTS");
	asserta(SIZE(Mega::m_Seqs) == QueryCount);

	const uint LA = M.m_ColCount;
	const uint ThreadCount = GetRequestedThreadCount();
	vector<float> Scores(QueryCount);
	vector<string> A3MSeqs(QueryCount);

	uint Counter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int qi = 0; qi < (int) QueryCount; ++qi)
		{
		LOCK();
		ProgressStep(Counter++, QueryCount, "Searching");
		UNLOCK();

		XDPMem &Mem = GetDPMem();
		const vector<vector<byte> > &Q = Mega::GetProfile((uint) qi);
		const string &Seq = Mega::m_Seqs[(uint) qi];
		string Path;
		uint Loi = 0;
		uint Loj = 0;
		float Score;
		if (Local)
			{
			uint Leni, Lenj;
			Score = SWFast_MASM(Mem, M, Q, Loi, Loj, Leni, Lenj, Path);
			}
		else
			Score = NWFast_MASM_MegaProf(Mem, M, Q, Loj, Path);

		Scores[(uint) qi] = Score;
		PathToA3M(Seq, LA, Loi, Loj, Path, A3MSeqs[(uint) qi]);
		}

	double Cutoff = 1e-3;
	if (M.m_HasCalibrate)
		Cutoff = M.m_CalibCutoff;
	if (optset_pvalue)
		{
		if (!M.m_HasCalibrate)
			Die("masm_search: -pvalue requires a calibrated MASM");
		Cutoff = opt(pvalue);
		if (Cutoff <= 0 || Cutoff > 1)
			Die("masm_search: -pvalue must be >0 and <=1");
		}

	vector<double> PValues(QueryCount, 1.0);
	vector<bool> Keep(QueryCount, true);
	if (M.m_HasCalibrate)
		{
		if (M.m_CalibSamples == 0)
			Die("masm_search: calibrate samples is 0");
		uint Pass = 0;
		for (uint i = 0; i < QueryCount; ++i)
			{
			double P = exp(M.m_CalibIntercept +
			  M.m_CalibSlope*double(Scores[i])) /
			  double(M.m_CalibSamples);
			if (P > 1)
				P = 1;
			PValues[i] = P;
			Keep[i] = (P <= Cutoff);
			if (Keep[i])
				++Pass;
			}
		ProgressLog("P-value cutoff %.4g  hits %u / %u\n",
		  Cutoff, Pass, QueryCount);
		}

	if (optset_output)
		{
		FILE *fTSV = CreateStdioFile(opt(output));
		const string &MasmLabel = M.m_Label;
		for (uint i = 0; i < QueryCount; ++i)
			{
			if (!Keep[i])
				continue;
			if (M.m_HasCalibrate)
				fprintf(fTSV, "%s\t%s\t%.3g\t%.4g\n",
				  Mega::GetLabel(i).c_str(), MasmLabel.c_str(),
				  Scores[i], PValues[i]);
			else
				fprintf(fTSV, "%s\t%s\t%.3g\n",
				  Mega::GetLabel(i).c_str(), MasmLabel.c_str(),
				  Scores[i]);
			}
		CloseStdioFile(fTSV);
		}

	if (optset_a3m)
		{
		FILE *fA3M = CreateStdioFile(opt(a3m));
		for (uint i = 0; i < QueryCount; ++i)
			{
			if (!Keep[i])
				continue;
			fprintf(fA3M, ">%s\n%s\n",
			  Mega::GetLabel(i).c_str(), A3MSeqs[i].c_str());
			}
		CloseStdioFile(fA3M);
		}
	}

void cmd_strumm_search()
	{
	if (optset_strumm)
		{
		optset_masm = true;
		opt_masm = opt_strumm;
		optused_masm = true;
		optused_strumm = true;
		}
	cmd_masm_search();
	}
