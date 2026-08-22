#pragma once

#include "myutils.h"
#include "msa_prep_colmap.h"
#include "msa.h"

enum class MSAPrepColMode
	{
	Auto,
	OccupancyOnly,
	LowerCaseInserts
	};

struct MSAPrepParams
	{
	MSAPrepColMode ColMode = MSAPrepColMode::Auto;
	double MaxGapFract = 0.5;
	uint MinCoreBlockCols = 8;
	double InsertMinorityFrac = 0.15;
	uint InsertDiscardLen = 30;
	uint InsertCapLen = 10;
	uint MergeGapCols = 3;
	uint SpacerMinCols = 2;
	uint SpacerMaxCols = 8;
	double SpacerMinOcc = 0.05;
	double SpacerMaxOcc = 0.35;
	};

struct MSAPrepBlock
	{
	uint BlockId = 0;
	uint OrigColLo = 0;
	uint OrigColHi = 0;
	uint CleanColLo = 0;
	uint CleanColHi = 0;
	};

struct MSAPrepInsertRegion
	{
	uint RegionId = 0;
	uint OrigColLo = 0;
	uint OrigColHi = 0;
	enum class Policy { Discard, Compress, Spacer } m_Policy = Policy::Discard;
	uint CleanWidth = 0;
	uint CleanColLo = UINT_MAX;
	vector<uint> OrigCols;
	};

struct MSAPrepResult
	{
	MSA CleanedMSA;
	vector<string> FullUngappedSeqs;
	vector<MSAPrepBlock> Blocks;
	vector<MSAPrepInsertRegion> InsertRegions;
	vector<MSAPrepColMapEntry> ColMap;
	uint OrigColCount = 0;
	uint CleanColCount = 0;
	};

MSAPrepResult MSAPrep(const MSA &Input, const MSAPrepParams &Params);

MSAPrepParams GetMSAPrepParamsFromOpts();

void WriteMSAPrepBlocksFile(const string &FileName,
  const MSAPrepResult &Result, const MSAPrepParams &Params);

void WriteMSAPrepMapFile(const string &FileName,
  const MSAPrepResult &Result);

void WriteMSAPrepBlocksJalView(const string &FileName,
  const MSA &OrigAln, const MSAPrepResult &Result);

void MSAToMultiSequence(const MSA &Aln, class MultiSequence &MS);
