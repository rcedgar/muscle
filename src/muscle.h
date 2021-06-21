#pragma once

#include "myutils.h"
#include "probcons.h"
#include "treeperm.h"

#if	DEBUG && !_DEBUG
#define _DEBUG	1
#endif

#if	_DEBUG && !DEBUG
#define DEBUG	1
#endif

#if	_MSC_VER
#define TIMING	0
#endif

#define VER_3_52	0

#ifdef	_MSC_VER	// Miscrosoft compiler
#pragma warning(disable : 4800)	// int-bool conversion
#pragma warning(disable : 4996)	// deprecated names like strdup, isatty.
#endif

#define SHORT_VERSION "5.0"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>

#define DOUBLE_AFFINE	0
#define SINGLE_AFFINE	1
#define PAF				0

#include "types.h"
#include "intmath.h"
#include "alpha.h"
#include "params.h"
#include "msa.h"

#ifndef _WIN32
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define	_snprintf snprintf
#define _fsopen(name, mode, share)	fopen((name), (mode))
#endif

typedef unsigned uint;

#if	DEBUG
#undef	assert
#define assert(b)	Call_MY_ASSERT(__FILE__, __LINE__, b, #b)
void Call_MY_ASSERT(const char *file, int line, bool b, const char *msg);
#else
#undef	assert
#define assert(exp)     ((void)0)
#endif

extern int g_argc;
extern char **g_argv;

#define Rotate(a, b, c)	{ SCORE *tmp = a; a = b; b = c; c = tmp; }

const double VERY_LARGE_DOUBLE = 1e20;

extern unsigned g_uTreeSplitNode1;
extern unsigned g_uTreeSplitNode2;

// Number of elements in array a[]
#define countof(a)	(sizeof(a)/sizeof(a[0]))

// Maximum of two of any type
#define	Max2(a, b)			((a) > (b) ? (a) : (b))

// Maximum of three of any type
#define	Max3(a, b, c)		Max2(Max2(a, b), c)

// Minimum of two of any type
#define Min2(a, b)		((a) < (b) ? (a) : (b))

// Maximum of four of any type
#define Max4(a, b, c, d)	Max2(Max2(a, b), Max2(c, d))

const double VERY_NEGATIVE_DOUBLE = -9e29;
const float VERY_NEGATIVE_FLOAT = (float) -9e29;

const double BLOSUM_DIST = 0.62;	// todo settable

// insane value for uninitialized variables
const unsigned uInsane = 8888888;
const int iInsane = 8888888;
const SCORE scoreInsane = 8888888;
const char cInsane = (char) 0xcd;		// int 3 instruction, used e.g. for unint. memory
const double dInsane = VERY_NEGATIVE_DOUBLE;
const float fInsane = VERY_NEGATIVE_FLOAT;
const char INVALID_STATE = '*';
const BASETYPE BTInsane = (BASETYPE) dInsane;
const WEIGHT wInsane = BTInsane;

extern double g_dNAN;

extern unsigned long g_tStart;

void Quit(const char szFormat[], ...);
void TrimBlanks(char szStr[]);
void TrimLeadingBlanks(char szStr[]);
void TrimTrailingBlanks(char szStr[]);
void Log(const char szFormat[], ...);
bool Verbose();
const char *ScoreToStr(SCORE Score);
const char *ScoreToStrL(SCORE Score);
SCORE StrToScore(const char *pszStr);
void Break();

double VecSum(const double v[], unsigned n);
bool IsValidInteger(const char *Str);
bool IsValidSignedInteger(const char *Str);
bool IsValidIdentifier(const char *Str);
bool IsValidFloatChar(char c);
bool isident(char c);
bool isidentf(char c);

void TreeFromSeqVect(const SeqVect &c, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root);
void TreeFromMSA(const MSA &msa, Tree &tree, CLUSTER Cluster,
  DISTANCE Distance, ROOT Root);

void StripGaps(char szStr[]);
void StripWhitespace(char szStr[]);
const char *GetTimeAsStr();
unsigned CalcBLOSUMWeights(MSA &Aln, ClusterTree &BlosumCluster);
void CalcGSCWeights(MSA &Aln, const ClusterTree &BlosumCluster);
void AssertNormalized(const PROB p[]);
void AssertNormalizedOrZero(const PROB p[]);
void AssertNormalized(const double p[]);
bool VectorIsZero(const double dValues[], unsigned n);
void VectorSet(double dValues[], unsigned n, double d);
bool VectorIsZero(const float dValues[], unsigned n);
void VectorSet(float dValues[], unsigned n, float d);

double pow2(double x);
double lnTolog2(double ln);

double lp2(double x);
SCORE SumLog(SCORE x, SCORE y);
SCORE SumLog(SCORE x, SCORE y, SCORE z);
SCORE SumLog(SCORE w, SCORE x, SCORE y, SCORE z);

double lp2Fast(double x);
double SumLogFast(double x, double y);
double SumLogFast(double x, double y, double z);
double SumLogFast(double w, double x, double y, double z);

void chkmem(const char szMsg[] = "");

void Normalize(PROB p[], unsigned n);
void Normalize(PROB p[], unsigned n, double dRequiredTotal);
void NormalizeUnlessZero(PROB p[], unsigned n);

void DebugPrintf(const char szFormat[], ...);
void SetListFileName(const char *ptrListFileName, bool bAppend);
void ModelFromAlign(const char *strInputFileName, const char *strModelFileName,
  double dMaxNIC);
double GetMemUseMB();
double GetRAMSizeMB();
double GetPeakMemUseMB();
void CheckMemUse();
const char *ElapsedTimeAsString();
char *SecsToHHMMSS(long lSecs, char szStr[]);
double GetCPUGHz();
SCORE GetBlosum62(unsigned uLetterA, unsigned uLetterB);
SCORE GetBlosum62d(unsigned uLetterA, unsigned uLetterB);
SCORE GetBlosum50(unsigned uLetterA, unsigned uLetterB);
void AssertNormalizedDist(const PROB p[], unsigned N);
void CmdLineError(const char *Format, ...);
void Fatal(const char *Format, ...);
void InitCmd();
void ExecCommandLine(int argc, char *argv[]);
void DoCmd();
void SetLogFile();
void NameFromPath(const char szPath[], char szName[], unsigned uBytes);
char *strsave(const char *s);
void DistKmer20_3(const SeqVect &v, DistFunc &DF);
void DistKbit20_3(const SeqVect &v, DistFunc &DF);
void DistKmer6_6(const SeqVect &v, DistFunc &DF);
void DistKmer4_6(const SeqVect &v, DistFunc &DF);
void DistPWKimura(const SeqVect &v, DistFunc &DF);
void FastDistKmer(const SeqVect &v, DistFunc &DF);
void DistUnaligned(const SeqVect &v, DISTANCE DistMethod, DistFunc &DF);
double PctIdToMAFFTDist(double dPctId);
double KimuraDist(double dPctId);
void SetFastParams();
void AssertProfsEq(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB);
void ValidateMuscleIds(const MSA &msa);
void ValidateMuscleIds(const Tree &tree);
void TraceBackToPath(int **TraceBack, unsigned uLengthA,
  unsigned uLengthB, PWPath &Path);
void BitTraceBack(char **TraceBack, unsigned uLengthA, unsigned uLengthB,
  char LastEdge, PWPath &Path);
SCORE AlignTwoMSAs(const MSA &msa1, const MSA &msa2, MSA &msaOut, PWPath &Path,
  bool bLockLeft = false, bool bLockRight = false);
SCORE AlignTwoProfs(
  const ProfPos *PA, unsigned uLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uLengthB, WEIGHT wB,
  PWPath &Path, ProfPos **ptrPout, unsigned *ptruLengthOut);
void AlignTwoProfsGivenPath(const PWPath &Path,
  const ProfPos *PA, unsigned uLengthA, WEIGHT wA,
  const ProfPos *PB, unsigned uLengthB, WEIGHT wB,
  ProfPos **ptrPOut, unsigned *ptruLengthOut);
void AlignTwoMSAsGivenPathSW(const PWPath &Path, const MSA &msaA, const MSA &msaB,
  MSA &msaCombined);
void AlignTwoMSAsGivenPath(const PWPath &Path, const MSA &msaA, const MSA &msaB,
  MSA &msaCombined);
SCORE FastScorePath2(const ProfPos *PA, unsigned uLengthA,
  const ProfPos *PB, unsigned uLengthB, const PWPath &Path);
SCORE GlobalAlignDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE GlobalAlignSimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE GlobalAlignSP(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE GlobalAlignSPN(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
SCORE GlobalAlignLE(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
void CalcThreeWayWeights(const Tree &tree, unsigned uNode1, unsigned uNode2,
  WEIGHT *Weights);
SCORE GlobalAlignSS(const Seq &seqA, const Seq &seqB, PWPath &Path);
bool RefineHoriz(MSA &msaIn, const Tree &tree, unsigned uIters, bool bLockLeft, bool bLockRight);
bool RefineVert(MSA &msaIn, const Tree &tree, unsigned uIters);
SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);

void SetInputFileName(const char *pstrFileName);
void SetIter(unsigned uIter);
void IncIter();
void SetMaxIters(unsigned uMaxIters);
void Progress(unsigned uStep, unsigned uTotalSteps);
void Progress(const char *szFormat, ...);
void SetStartTime();
void ProgressStepsDone();
void SetProgressDesc(const char szDesc[]);
void ShowSeqStats(unsigned uSeqCount, unsigned uMinL, unsigned uMaxL, unsigned uAvgL);

void SetNewHandler();
void SaveCurrentAlignment();
void SetCurrentAlignment(MSA &msa);
void SetOutputFileName(const char *out);

#if	DEBUG
void SetMuscleSeqVect(SeqVect &v);
void SetMuscleInputMSA(MSA &msa);
void ValidateMuscleIds(const MSA &msa);
void ValidateMuscleIds(const Tree &tree);
#else
#define SetMuscleSeqVect(x)		/* empty */
#define SetMuscleInputMSA(x)	/* empty */
#define ValidateMuscleIds(x)	/* empty */
#endif

void SetParams();

void SortCounts(const FCOUNT fcCounts[], unsigned SortOrder[]);
unsigned ResidueGroupFromFCounts(const FCOUNT fcCounts[]);
FCOUNT SumCounts(const FCOUNT Counts[]);

bool FlagOpt(const char *Name);
const char *ValueOpt(const char *Name);
void DoMuscle();
void ProfDB();
void DoSP();
void ProgAlignSubFams();
void Run();
void LogParams();
void OnException();
void SetSeqWeightMethod(SEQWEIGHT Method);
SEQWEIGHT GetSeqWeightMethod();
WEIGHT GetMuscleSeqWeightById(unsigned uId);
void ListDiagSavings();
void CheckMaxTime();
const char *MaxSecsToStr();
unsigned long GetStartTime();

void ProgressiveAlign(const SeqVect &v, const Tree &GuideTree, MSA &a);

void MakeRootMSA(const SeqVect &v, const Tree &GuideTree, ProgNode Nodes[], MSA &a);
void MakeRootMSABrenner(SeqVect &v, const Tree &GuideTree, ProgNode Nodes[], MSA &a);

void Refine();
void Local();
void Profile();
void PPScore();
void UPGMA2(const DistCalc &DC, Tree &tree, LINKAGE Linkage);

char *GetFastaSeq(FILE *f, unsigned *ptrSeqLength, char **ptrLabel,
  bool DeleteGaps = true);
SCORE SW(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
void TraceBackSW(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  unsigned uPrefixLengthAMax, unsigned uPrefixLengthBMax, PWPath &Path);
void DiffPaths(const PWPath &p1, const PWPath &p2, unsigned Edges1[],
  unsigned *ptruDiffCount1, unsigned Edges2[], unsigned *ptruDiffCount2);
void SetPPScore(bool bRespectFlagOpts = true);
void SetPPScore(PPSCORE p);
SCORE GlobalAlignDimer(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
bool MissingCommand();
void Credits();
void ProfileProfile(MSA &msa1, MSA &msa2, MSA &msaOut);
void MHackStart(SeqVect &v);
void MHackEnd(MSA &msa);
void WriteScoreFile(const MSA &msa);
char ConsensusChar(const ProfPos &PP);
void Stabilize(const MSA &msa, MSA &msaStable);
void MuscleOutput(MSA &msa);
PTR_SCOREMATRIX ReadMx(TextFile &File);
void MemPlus(size_t Bytes, char *Where);
void MemMinus(size_t Bytes, char *Where);

void SetAlphaFromSeqVect(SeqVect &SV);
void ShowSeqStats(const SeqVect &SV);
void GetConsensusSequence(const MultiSequence &MSA, string &Seq);
void TreeFromSeqVect_EA(const SeqVect &SV, Tree &SplitGuideTree2);
Sequence *SeqToSequence(const Seq &inseq, uint Index);
void GetLabelToIndex(const SeqVect &SV, map<string, uint> &LabelToIndex);
void GetSeqsByLabels(const SeqVect &SV, const vector<string> &Labels,
  const map<string, uint> &LabelToIndex, MultiSequence &Seqs);
void ProgressLogInputSummary(const string &FileName, const MultiSequence &Seqs);
void RunMTProbcons(MultiSequence &InputSeqs);
void RunSuper4(MultiSequence &InputSeqs, MultiSequence &MSA, TREEPERM TreePerm);

void CalcEADistMx(FILE *f, MultiSequence* sequences,
  vector<vector<float> > &DistMx);
void PermuteTree(const Tree &InputTree,
  Tree &TreeABC, Tree &TreeACB, Tree &TreeBCA,
  vector<string> &LabelsA, vector<string> &LabelsB, vector<string> &LabelsC);
void PermTree(Tree &InputTree, TREEPERM TP);
void StringsToFile(const string &FileName, const vector<string> &v);
void MakeFileName(const string &Pattern, uint N, string &FileName);
MultiSequence &LoadGlobalInputMS(const string &FileName);
const MultiSequence &GetGlobalInputMS();
uint GetGlobalMSSeqCount();
uint GetGSICount();
uint GetAssertSameSeqsOkCount();
const Sequence &GetGlobalInputSeq(uint GSI);
const string &GetGlobalInputSeqLabel(uint GSI);
void CharVecToStr(const vector<char> &Vec, string &Str);
void LogAln(const Sequence &X, const Sequence &Y, const string &PathXY);
void ReadStringsFromFile(const string &FileName,
  vector<string> &Strings);
void GetGuideTreeJoinOrder(const Tree &GuideTree,
  const map<string, uint> &LabelToIndex,
  vector<uint> &Indexes1, vector<uint> &Indexes2);
void ValidateJoinOrder(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2);

extern FILE *g_fColEA;
extern MultiSequence *g_RefMSA;

void CalcEA_MSA(const MultiSequence &MSA, 
  vector<vector<SparseMatrix*> > &sparseMatrices,
  vector<float> &ColToEA);
void WriteColEA(FILE *f, const MultiSequence &MSA,
  vector<float> &ColToEA);
float RefEA(const MultiSequence &TestMSA, const MultiSequence &RefMSA,
  const vector<float > &ColToEA);

void _AssertSameLabels(const char *File, uint Line,
  MultiSequence &MS);
void _AssertSameSeqs(const char *File, uint Line, 
  MultiSequence &MS1, MultiSequence &MS2);
void _AssertSameSeqsVec(const char *File, uint Line, 
  MultiSequence &MS, vector<MultiSequence *> &v);
void _AssertSameSeqsJoin(const char *File, uint Line, 
  MultiSequence &MS1, MultiSequence &MS2, MultiSequence &MS12);
void _AssertSeqsEqInput(const char *File, uint Line,
  const MultiSequence &MS);

#define AssertSeqsEqInput(MS)		_AssertSeqsEqInput(__FILE__, __LINE__, MS)
#define AssertSameLabels(MS)		_AssertSameLabels(__FILE__, __LINE__, MS)
#define AssertSameSeqs(MS1, MS2)	_AssertSameSeqs(__FILE__, __LINE__, MS1, MS2)
#define AssertSameSeqsVec(MS, v)	_AssertSameSeqsVec(__FILE__, __LINE__, MS, v)
#define AssertSameSeqsJoin(MS1, MS2, MS12)	_AssertSameSeqsJoin(__FILE__, __LINE__, MS1, MS2, MS12)
