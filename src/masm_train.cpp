#include "muscle.h"
#include "msa.h"
#include "msa_prep.h"
#include "masm.h"
#include "mega.h"
#include "multisequence.h"

void cmd_masm_stats()
	{
	MASM M;
	M.FromFile(g_Arg1);
	ProgressLog("%10u  Sequences\n", M.m_SeqCount);
	ProgressLog("%10u  Columns\n", M.m_ColCount);
	ProgressLog("%10u  Features ", M.m_FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < M.m_FeatureCount; ++ FeatureIdx)
		ProgressLog(" %s/%u",
		  M.m_FeatureNames[FeatureIdx].c_str(), 
		  M.m_AlphaSizes[FeatureIdx]); 
	ProgressLog("\n");
	if (M.m_SeedMSA.GetSeqCount() > 0)
		ProgressLog("%10u  Seed MSA seqs (%u cols)\n",
		  M.m_SeedMSA.GetSeqCount(), M.m_SeedMSA.GetColCount());
	}

void cmd_masm_train()
	{
	const string &AlnFN = g_Arg1;
	const string &StructsFN = opt(input);

	Mega::RejectLegacyMega(StructsFN);
	Mega::FromStructs(StructsFN);

	MSA OrigMSA;
	OrigMSA.FromFASTAFile_PreserveCase(AlnFN);

	string Label;
	if (optset_label)
		Label = opt(label);
	else
		Label = string(BaseName(AlnFN.c_str()));

	float GapOpen = 4;
	float GapExt = 0.5;
	if (optset_gapopen)
		GapOpen = (float) opt(gapopen);
	if (optset_gapext)
		GapExt = (float) opt(gapext);

	const MSAPrepResult *PrepPtr = 0;
	MSAPrepResult Prep;
	if (!opt(noprep))
		{
		MSAPrepParams Params = GetMSAPrepParamsFromOpts();
		Prep = MSAPrep(OrigMSA, Params);
		PrepPtr = &Prep;
		}

	if (optset_jalview_features)
		{
		if (PrepPtr == 0)
			Die("-jalview_features requires prep (omit -noprep)");
		WriteMSAPrepBlocksJalView(opt(jalview_features), OrigMSA, *PrepPtr);
		}

	string MapFN;
	if (optset_map)
		{
		if (PrepPtr == 0)
			Die("-map requires prep (omit -noprep)");
		MapFN = opt(map);
		}
	else if (PrepPtr != 0 && optset_output)
		MapFN = string(opt(output)) + ".map";
	if (!MapFN.empty())
		WriteMSAPrepMapFile(MapFN, *PrepPtr);

	const MSA *TrainMSA = (PrepPtr != 0 ? &PrepPtr->CleanedMSA : &OrigMSA);
	const vector<MSAPrepColMapEntry> *ColMap =
	  (PrepPtr != 0 ? &PrepPtr->ColMap : 0);
	const vector<string> *FullUngapped =
	  (PrepPtr != 0 ? &PrepPtr->FullUngappedSeqs : 0);

	MultiSequence TrainAln;
	MSAToMultiSequence(*TrainMSA, TrainAln);
	if (optset_seedmsaout)
		TrainAln.ToFasta(opt(seedmsaout));

	MASM M;
	M.FromMSA(TrainAln, Label, GapOpen, GapExt, ColMap, FullUngapped);

	if (!opt(nocalibrate))
		{
		if (opt(local) && opt(global))
			Die("masm_train: specify at most one of -local or -global");
		const bool Local = opt(local);

		bool Denovo = opt(denovo);
		bool Decoy = optset_decoy;
		bool Shatter = opt(shatter);
		uint FPModeCount = 0;
		if (Denovo)
			++FPModeCount;
		if (Decoy)
			++FPModeCount;
		if (Shatter)
			++FPModeCount;
		if (FPModeCount > 1)
			Die("masm_train: specify at most one of -denovo, -decoy or -shatter");
		if (FPModeCount == 0)
			{
			if (!MapFN.empty())
				Denovo = true;
			else
				Shatter = true;
			}
		if (Denovo && MapFN.empty())
			Die("masm_train -denovo: map file required (omit -noprep or set -map)");
		if (Decoy && opt(decoy)[0] == 0)
			Die("masm_train: -decoy requires a MASM file");

		const uint N = optset_n ? opt(n) : 10000;
		uint ShatterMin = optset_shatter_min ? opt(shatter_min) : 5;
		uint ShatterMax = optset_shatter_max ? opt(shatter_max) : 15;
		string HistTSV;
		if (optset_tsvout)
			HistTSV = opt(tsvout);

		CalibrateMASM(M, Local, N,
		  Denovo, MapFN,
		  Decoy, Decoy ? opt(decoy) : "",
		  Shatter, ShatterMin, ShatterMax,
		  HistTSV);
		}

	M.ToFile(opt(output));
	ProgressLog("Wrote MASM %u cols to %s\n",
	  M.GetColCount(), opt(output));
	}

void cmd_strumm_build()
	{
	cmd_masm_train();
	}