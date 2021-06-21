#include "muscle.h"
#include "objscore.h"
#include "profile.h"
#include "enumopts.h"

const double DEFAULT_MAX_MB_FRACT = 0.8;

SCORE g_scoreCenter = 0;
SCORE g_scoreGapExtend = 0;
SCORE g_scoreGapOpen2 = MINUS_INFINITY;
SCORE g_scoreGapExtend2 = MINUS_INFINITY;
SCORE g_scoreGapAmbig = 0;
SCORE g_scoreAmbigFactor = 0;

extern SCOREMATRIX VTML_LA;
extern SCOREMATRIX PAM200;
extern SCOREMATRIX PAM200NoCenter;
extern SCOREMATRIX VTML_SP;
extern SCOREMATRIX VTML_SPNoCenter;
extern SCOREMATRIX NUC_SP;

PTR_SCOREMATRIX g_ptrScoreMatrix;

const char *g_pstrInFileName = "-";
const char *g_pstrOutFileName = "-";
const char *g_pstrFASTAOutFileName = 0;
const char *g_pstrMSFOutFileName = 0;
const char *g_pstrClwOutFileName = 0;
const char *g_pstrClwStrictOutFileName = 0;
const char *g_pstrHTMLOutFileName = 0;
const char *g_pstrPHYIOutFileName = 0;
const char *g_pstrPHYSOutFileName = 0;

const char *g_pstrFileName1 = 0;
const char *g_pstrFileName2 = 0;

const char *g_pstrSPFileName = 0;
const char *g_pstrMatrixFileName = 0;

const char *g_pstrUseTreeFileName = 0;
bool g_bUseTreeNoWarn = false;

const char *g_pstrComputeWeightsFileName;
const char *g_pstrScoreFileName;

const char *g_pstrProf1FileName = 0;
const char *g_pstrProf2FileName = 0;

unsigned g_uSmoothWindowLength = 7;
unsigned g_uAnchorSpacing = 32;
unsigned g_uMaxTreeRefineIters = 1;

unsigned g_uRefineWindow = 200;
unsigned g_uWindowFrom = 0;
unsigned g_uWindowTo = 0;
unsigned g_uSaveWindow = uInsane;
unsigned g_uWindowOffset = 0;

unsigned g_uMaxSubFamCount = 5;

unsigned g_uHydrophobicRunLength = 5;
float g_dHydroFactor = (float) 1.2;

unsigned g_uMinDiagLength = 24;	// TODO alpha -- should depend on alphabet?
unsigned g_uMaxDiagBreak = 1;
unsigned g_uDiagMargin = 5;

float g_dSUEFF = (float) 0.1;

bool g_bPrecompiledCenter = true;
bool g_bNormalizeCounts = false;
bool g_bDiags1 = false;
bool g_bDiags2 = false;
bool g_bAnchors = true;
bool g_bQuiet = false;
bool g_bVerbose = false;
bool g_bRefine = false;
bool g_bRefineW = false;
bool g_bProfDB = false;
bool g_bLow = false;
bool g_bSW = false;
bool g_bClusterOnly = false;
bool g_bProfile = false;
bool g_bPPScore = false;
bool g_bBrenner = false;
bool g_bDimer = false;
bool g_bVersion = false;
bool g_bStable = false;
bool g_bFASTA = false;
bool g_bTomHydro = false;
bool g_bMakeTree = false;
bool g_bMSF = false;
bool g_bAln = false;
bool g_bClwStrict = false;
bool g_bHTML = false;
bool g_bPHYI = false;
bool g_bPHYS = false;

unsigned g_uMaxIters = 8;
unsigned long g_ulMaxSecs = 0;
unsigned g_uMaxMB = 500;

PPSCORE g_PPScore = PPSCORE_LE;
OBJSCORE g_ObjScore = OBJSCORE_SPM;

SEQWEIGHT g_SeqWeight1 = SEQWEIGHT_ClustalW;
SEQWEIGHT g_SeqWeight2 = SEQWEIGHT_ClustalW;

DISTANCE g_Distance1 = DISTANCE_Kmer6_6;
DISTANCE g_Distance2 = DISTANCE_PctIdKimura;

CLUSTER g_Cluster1 = CLUSTER_UPGMB;
CLUSTER g_Cluster2 = CLUSTER_UPGMB;

ROOT g_Root1 = ROOT_FromClustering;
ROOT g_Root2 = ROOT_FromClustering;

bool g_bDiags;

SEQTYPE g_SeqType = SEQTYPE_Auto;

TERMGAPS g_TermGaps = TERMGAPS_Half;

//------------------------------------------------------
// These parameters depend on the chosen prof-prof
// score (g_PPScore), initialized to "Undefined".
float g_dSmoothScoreCeil = fInsane;
float g_dMinBestColScore = fInsane;
float g_dMinSmoothScore = fInsane;
SCORE g_scoreGapOpen = fInsane;
//------------------------------------------------------

const char *MaxSecsToStr()
	{
	if (0 == g_ulMaxSecs)
		return "(No limit)";
	return SecsToStr(g_ulMaxSecs);
	}

void LogParams()
	{
	Log("Profile-profile score    %s\n", PPSCOREToStr(g_PPScore));
	Log("Max iterations           %u\n", g_uMaxIters);
	Log("Max trees                %u\n", g_uMaxTreeRefineIters);
	Log("Max time                 %s\n", MaxSecsToStr());
	Log("Max MB                   %u\n", g_uMaxMB);
	Log("Gap open                 %g\n", g_scoreGapOpen);
	Log("Gap extend (dimer)       %g\n", g_scoreGapExtend);
	Log("Gap ambig factor         %g\n", g_scoreAmbigFactor);
	Log("Gap ambig penalty        %g\n", g_scoreGapAmbig);
	Log("Center (LE)              %g\n", g_scoreCenter);
	Log("Term gaps                %s\n", TERMGAPSToStr(g_TermGaps));

	Log("Smooth window length     %u\n", g_uSmoothWindowLength);
	Log("Refine window length     %u\n", g_uRefineWindow);
	Log("Min anchor spacing       %u\n", g_uAnchorSpacing);
	Log("Min diag length (lambda) %u\n", g_uMinDiagLength);
	Log("Diag margin (mu)         %u\n", g_uDiagMargin);
	Log("Min diag break           %u\n", g_uMaxDiagBreak);
	Log("Hydrophobic window       %u\n", g_uHydrophobicRunLength);

	Log("Hydrophobic gap factor   %g\n", g_dHydroFactor);
	Log("Smooth score ceiling     %g\n", g_dSmoothScoreCeil);
	Log("Min best col score       %g\n", g_dMinBestColScore);
	Log("Min anchor score         %g\n", g_dMinSmoothScore);
	Log("SUEFF                    %g\n", g_dSUEFF);

	Log("Brenner root MSA         %s\n", BoolToStr(g_bBrenner));
	Log("Normalize counts         %s\n", BoolToStr(g_bNormalizeCounts));
	Log("Diagonals (1)            %s\n", BoolToStr(g_bDiags1));
	Log("Diagonals (2)            %s\n", BoolToStr(g_bDiags2));
	Log("Anchors                  %s\n", BoolToStr(g_bAnchors));
	Log("MSF output format        %s\n", BoolToStr(g_bMSF));
	Log("Phylip interleaved       %s\n", BoolToStr(g_bPHYI));
	Log("Phylip sequential        %s\n", BoolToStr(g_bPHYS));
	Log("ClustalW output format   %s\n", BoolToStr(g_bAln));
	Log("Quiet                    %s\n", BoolToStr(g_bQuiet));
	Log("Refine                   %s\n", BoolToStr(g_bRefine));
	Log("ProdfDB                  %s\n", BoolToStr(g_bProfDB));
	Log("Low complexity profiles  %s\n", BoolToStr(g_bLow));

	Log("Objective score          %s\n", OBJSCOREToStr(g_ObjScore));

	Log("Distance method (1)      %s\n", DISTANCEToStr(g_Distance1));
	Log("Clustering method (1)    %s\n", CLUSTERToStr(g_Cluster1));
	Log("Root method (1)          %s\n", ROOTToStr(g_Root1));
	Log("Sequence weighting (1)   %s\n", SEQWEIGHTToStr(g_SeqWeight1));

	Log("Distance method (2)      %s\n", DISTANCEToStr(g_Distance2));
	Log("Clustering method (2)    %s\n", CLUSTERToStr(g_Cluster2));
	Log("Root method (2)          %s\n", ROOTToStr(g_Root2));
	Log("Sequence weighting (2)   %s\n", SEQWEIGHTToStr(g_SeqWeight2));

	Log("\n");
	}

static void SetDefaultsLE()
	{
	g_ptrScoreMatrix = &VTML_LA;

	g_scoreGapOpen = (SCORE) -2.9;
	g_scoreCenter = (SCORE) -0.52;

	g_bNormalizeCounts = true;

	g_dSmoothScoreCeil = 3.0;
	g_dMinBestColScore = 2.0;
	g_dMinSmoothScore = 1.0;

	g_Distance1 = DISTANCE_Kmer6_6;
	g_Distance2 = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSP()
	{
	g_ptrScoreMatrix = &PAM200;

	g_scoreGapOpen = -1439;
	g_scoreCenter = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 200.0;
	g_dMinBestColScore = 300.0;
	g_dMinSmoothScore = 125.0;

	g_Distance1 = DISTANCE_Kmer6_6;
	g_Distance2 = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSV()
	{
	g_ptrScoreMatrix = &VTML_SP;

	g_scoreGapOpen = -300;
	g_scoreCenter = 0.0;	// center pre-added into score mx

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 90.0;
	g_dMinBestColScore = 130.0;
	g_dMinSmoothScore = 40.0;

	g_Distance1 = DISTANCE_Kmer6_6;
	g_Distance2 = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSPN_DNA()
	{
	g_ptrScoreMatrix = &NUC_SP;

	g_scoreGapOpen = -400;
	g_scoreCenter = 0.0;	// center pre-added into score mx
	g_scoreGapExtend = 0.0;

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 999.0;		// disable
	g_dMinBestColScore = 90;
	g_dMinSmoothScore = 90;

	g_Distance1 = DISTANCE_Kmer4_6;
	g_Distance2 = DISTANCE_PctIdKimura;
	}

static void SetDefaultsSPN_RNA()
	{
	g_ptrScoreMatrix = &NUC_SP;

	g_scoreGapOpen = -420;
	g_scoreCenter = -300;	// total center = NUC_EXTEND - 300 
	g_scoreGapExtend = 0.0;

	g_bNormalizeCounts = false;

	g_dSmoothScoreCeil = 999.0;		// disable
	g_dMinBestColScore = 90;
	g_dMinSmoothScore = 90;

	g_Distance1 = DISTANCE_Kmer4_6;
	g_Distance2 = DISTANCE_PctIdKimura;
	}

static void FlagParam(const char *OptName, bool *ptrParam, bool bValueIfFlagSet)
	{
	bool bIsSet = FlagOpt(OptName);
	if (bIsSet)
		*ptrParam = bValueIfFlagSet;
	}

static void StrParam(const char *OptName, const char **ptrptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrptrParam = opt;
	}

static void FloatParam(const char *OptName, float *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = (float) atof(opt);
	}

static void UintParam(const char *OptName, unsigned *ptrParam)
	{
	const char *opt = ValueOpt(OptName);
	if (0 != opt)
		*ptrParam = atou(opt);
	}

static void EnumParam(const char *OptName, EnumOpt *Opts, int *Param)
	{
	const char *Value = ValueOpt(OptName);
	if (0 == Value)
		return;

	for (;;)
		{
		if (0 == Opts->pstrOpt)
			Quit("Invalid parameter -%s %s", OptName, Value);
		if (0 == stricmp(Value, Opts->pstrOpt))
			{
			*Param = Opts->iValue;
			return;
			}
		++Opts;
		}
	}

static void SetPPDefaultParams()
	{
	switch (g_PPScore)
		{
	case PPSCORE_SP:
		SetDefaultsSP();
		break;

	case PPSCORE_LE:
		SetDefaultsLE();
		break;

	case PPSCORE_SV:
		SetDefaultsSV();
		break;

	case PPSCORE_SPN:
		switch (g_Alpha)
			{
		case ALPHA_DNA:
			SetDefaultsSPN_DNA();
			break;
		case ALPHA_RNA:
			SetDefaultsSPN_RNA();
			break;
		default:
			Quit("Invalid alpha %d", g_Alpha);
			}
		break;

	default:
		Quit("Invalid g_PPScore");
		}
	}

static void SetPPCommandLineParams()
	{
	FloatParam("GapOpen", &g_scoreGapOpen);
	FloatParam("GapOpen2", &g_scoreGapOpen2);
	FloatParam("GapExtend", &g_scoreGapExtend);
	FloatParam("GapExtend2", &g_scoreGapExtend2);
	FloatParam("GapAmbig", &g_scoreAmbigFactor);
	FloatParam("Center", &g_scoreCenter);
	FloatParam("SmoothScoreCeil", &g_dSmoothScoreCeil);
	FloatParam("MinBestColScore", &g_dMinBestColScore);
	FloatParam("MinSmoothScore", &g_dMinSmoothScore);

	EnumParam("Distance", DISTANCE_Opts, (int *) &g_Distance1);
	EnumParam("Distance", DISTANCE_Opts, (int *) &g_Distance2);

	EnumParam("Distance1", DISTANCE_Opts, (int *) &g_Distance1);
	EnumParam("Distance2", DISTANCE_Opts, (int *) &g_Distance2);
	}

void SetPPScore(bool bRespectFlagOpts)
	{
	if (bRespectFlagOpts)
		{
		if (FlagOpt("SP"))
			g_PPScore = PPSCORE_SP;
		else if (FlagOpt("LE"))
			g_PPScore = PPSCORE_LE;
		else if (FlagOpt("SV"))
			g_PPScore = PPSCORE_SV;
		else if (FlagOpt("SPN"))
			g_PPScore = PPSCORE_SPN;
		}

	switch (g_PPScore)
		{
	case PPSCORE_LE:
	case PPSCORE_SP:
	case PPSCORE_SV:
		if (ALPHA_RNA == g_Alpha || ALPHA_DNA == g_Alpha)
			g_PPScore = PPSCORE_SPN;
		break;
	case PPSCORE_SPN:
		if (ALPHA_Amino == g_Alpha)
			g_PPScore = PPSCORE_LE;
		break;
		}

	SetPPDefaultParams();
	SetPPCommandLineParams();

	if (g_bVerbose)
		LogParams();
	}

void SetPPScore(PPSCORE p)
	{
	g_PPScore = p;
	SetPPScore(true);
	}

static void SetMaxSecs()
	{
	float fMaxHours = 0.0;
	FloatParam("MaxHours", &fMaxHours);
	if (0.0 == fMaxHours)
		return;
	g_ulMaxSecs = (unsigned long) (fMaxHours*60*60);
	}

static bool CanDoLowComplexity()
	{
	if (g_SeqWeight1 != SEQWEIGHT_ClustalW)
		return false;
	if (1 == g_uMaxIters)
		return true;
	return g_SeqWeight2 == SEQWEIGHT_ClustalW;
	}

bool MissingCommand()
	{
	if (strcmp(g_pstrInFileName, "-"))
		return false;
	if (0 != g_pstrFileName1)
		return false;
	if (0 != g_pstrSPFileName)
		return false;
	return true;
	}

void SetParams()
	{
	SetMaxSecs();

	StrParam("in", &g_pstrInFileName);
	StrParam("out", &g_pstrOutFileName);

	StrParam("FASTAOut", &g_pstrFASTAOutFileName);
	StrParam("ClwOut", &g_pstrClwOutFileName);
	StrParam("ClwStrictOut", &g_pstrClwStrictOutFileName);
	StrParam("HTMLOut", &g_pstrHTMLOutFileName);
	StrParam("PHYIOut", &g_pstrPHYIOutFileName);
	StrParam("PHYSOut", &g_pstrPHYSOutFileName);
	StrParam("MSFOut", &g_pstrMSFOutFileName);

	StrParam("in1", &g_pstrFileName1);
	StrParam("in2", &g_pstrFileName2);

	StrParam("Matrix", &g_pstrMatrixFileName);
	StrParam("SPScore", &g_pstrSPFileName);

	StrParam("UseTree_NoWarn", &g_pstrUseTreeFileName);
	if (0 != g_pstrUseTreeFileName)
		g_bUseTreeNoWarn = true;

	StrParam("UseTree", &g_pstrUseTreeFileName);
	StrParam("ComputeWeights", &g_pstrComputeWeightsFileName);
	StrParam("ScoreFile", &g_pstrScoreFileName);

	FlagParam("Diags1", &g_bDiags1, true);
	FlagParam("Diags2", &g_bDiags2, true);

	bool Diags = false;
	FlagParam("Diags", &Diags, true);
	if (Diags)
		{
		g_bDiags1 = true;
		g_bDiags2 = true;
		}

	FlagParam("Anchors", &g_bAnchors, true);
	FlagParam("NoAnchors", &g_bAnchors, false);

	FlagParam("Quiet", &g_bQuiet, true);
	FlagParam("Verbose", &g_bVerbose, true);
	FlagParam("Version", &g_bVersion, true);
	FlagParam("Stable", &g_bStable, true);
	FlagParam("Group", &g_bStable, false);
	FlagParam("Refine", &g_bRefine, true);
	FlagParam("RefineW", &g_bRefineW, true);
	FlagParam("ProfDB", &g_bProfDB, true);
	FlagParam("SW", &g_bSW, true);
	FlagParam("ClusterOnly", &g_bClusterOnly, true);
	FlagParam("Profile", &g_bProfile, true);
	FlagParam("PPScore", &g_bPPScore, true);
	FlagParam("Brenner", &g_bBrenner, true);
	FlagParam("Dimer", &g_bDimer, true);

	FlagParam("MSF", &g_bMSF, true);
	FlagParam("PHYI", &g_bPHYI, true);
	FlagParam("PHYS", &g_bPHYS, true);
	FlagParam("clw", &g_bAln, true);
	FlagParam("HTML", &g_bHTML, true);
	FlagParam("FASTA", &g_bFASTA, true);
	FlagParam("MakeTree", &g_bMakeTree, true);

	if (g_bStable)
		Quit("-stable not supported in this version of muscle");

	bool b = false;
	FlagParam("clwstrict", &b, true);
	if (b)
		{
		g_bAln = true;
		g_bClwStrict = true;
		}

	UintParam("MaxIters", &g_uMaxIters);
	UintParam("MaxTrees", &g_uMaxTreeRefineIters);
	UintParam("SmoothWindow", &g_uSmoothWindowLength);
	UintParam("RefineWindow", &g_uRefineWindow);
	UintParam("FromWindow", &g_uWindowFrom);
	UintParam("ToWindow", &g_uWindowTo);
	UintParam("SaveWindow", &g_uSaveWindow);
	UintParam("WindowOffset", &g_uWindowOffset);
	UintParam("AnchorSpacing", &g_uAnchorSpacing);
	UintParam("DiagLength", &g_uMinDiagLength);
	UintParam("DiagMargin", &g_uDiagMargin);
	UintParam("DiagBreak", &g_uMaxDiagBreak);
	UintParam("MaxSubFam", &g_uMaxSubFamCount);

	UintParam("Hydro", &g_uHydrophobicRunLength);
	FlagParam("TomHydro", &g_bTomHydro, true);
	if (g_bTomHydro)
		g_uHydrophobicRunLength = 0;

	FloatParam("SUEFF", &g_dSUEFF);
	FloatParam("HydroFactor", &g_dHydroFactor);

	EnumParam("ObjScore", OBJSCORE_Opts, (int *) &g_ObjScore);
	EnumParam("TermGaps", TERMGAPS_Opts, (int *) &g_TermGaps);

	EnumParam("Weight", SEQWEIGHT_Opts, (int *) &g_SeqWeight1);
	EnumParam("Weight", SEQWEIGHT_Opts, (int *) &g_SeqWeight2);

	EnumParam("Weight1", SEQWEIGHT_Opts, (int *) &g_SeqWeight1);
	EnumParam("Weight2", SEQWEIGHT_Opts, (int *) &g_SeqWeight2);

	EnumParam("Cluster", CLUSTER_Opts, (int *) &g_Cluster1);
	EnumParam("Cluster", CLUSTER_Opts, (int *) &g_Cluster2);

	EnumParam("Cluster1", CLUSTER_Opts, (int *) &g_Cluster1);
	EnumParam("Cluster2", CLUSTER_Opts, (int *) &g_Cluster2);

	EnumParam("Root1", ROOT_Opts, (int *) &g_Root1);
	EnumParam("Root2", ROOT_Opts, (int *) &g_Root2);

	EnumParam("SeqType", SEQTYPE_Opts, (int *) &g_SeqType);

	g_scoreGapAmbig = g_scoreGapOpen*g_scoreAmbigFactor;
	g_bLow = CanDoLowComplexity();

	if (g_bDimer)
		g_bPrecompiledCenter = false;

	UintParam("MaxMB", &g_uMaxMB);
	if (0 == ValueOpt("MaxMB"))
		g_uMaxMB = (unsigned) (GetRAMSizeMB()*DEFAULT_MAX_MB_FRACT);
	}
