#include "muscle.h"
#include "msa.h"
#include "msa_prep.h"
#include "masm.h"
#include "mega.h"
#include "multisequence.h"

static void TrainAndWriteMASM(const MSA &OrigMSA, const MSAPrepResult *Prep,
  const string &Label, float GapOpen, float GapExt, const string &OutFN)
	{
	const MSA *TrainMSA = &OrigMSA;
	const vector<MSAPrepColMapEntry> *ColMap = 0;
	const vector<string> *FullUngapped = 0;
	if (Prep != 0)
		{
		TrainMSA = &Prep->CleanedMSA;
		ColMap = &Prep->ColMap;
		FullUngapped = &Prep->FullUngappedSeqs;
		}

	MultiSequence TrainAln;
	MSAToMultiSequence(*TrainMSA, TrainAln);

	MASM M;
	M.FromMSA(TrainAln, Label, GapOpen, GapExt, ColMap, FullUngapped);
	M.ToFile(OutFN);
	ProgressLog("Wrote MASM %u cols to %s\n",
	  M.GetColCount(), OutFN.c_str());
	}

void cmd_msa_prep()
	{
	const string &AlnFN = g_Arg1;

	MSA InputMSA;
	InputMSA.FromFASTAFile_PreserveCase(AlnFN);

	MSAPrepParams Params = GetMSAPrepParamsFromOpts();
	MSAPrepResult Result = MSAPrep(InputMSA, Params);

	if (optset_output)
		Result.CleanedMSA.ToFASTAFile(opt(output));

	if (optset_blocks)
		WriteMSAPrepBlocksFile(opt(blocks), Result, Params);

	if (optset_map)
		WriteMSAPrepMapFile(opt(map), Result);

	if (optset_jalview_features)
		WriteMSAPrepBlocksJalView(opt(jalview_features), InputMSA, Result);

	if (optset_output_masm)
		{
		if (!optset_input)
			Die("-output_masm requires -input STRUCTS file");
		const string &StructsFN = opt(input);
		Mega::RejectLegacyMega(StructsFN);
		Mega::FromStructs(StructsFN);

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

		TrainAndWriteMASM(InputMSA, &Result, Label, GapOpen, GapExt,
		  opt(output_masm));
		}

	if (!optset_output && !optset_blocks && !optset_map &&
	  !optset_jalview_features && !optset_output_masm)
		Die("msa_prep: specify at least one of -output, -blocks, "
		  "-map, -jalview_features, -output_masm");
	}
