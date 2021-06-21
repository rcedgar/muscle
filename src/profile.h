#ifndef FastProf2_h
#define FastProf2_h

#include "msa.h"
#include "pwpath.h"
#include <math.h>	// for log function

class DiagList;
class WeightList;

struct ProfPos
	{
	bool m_bAllGaps;
	unsigned m_uSortOrder[21];
	FCOUNT m_fcCounts[20];
	FCOUNT m_LL;
	FCOUNT m_LG;
	FCOUNT m_GL;
	FCOUNT m_GG;
	SCORE m_AAScores[20];
	unsigned m_uResidueGroup;
	FCOUNT m_fOcc;
	FCOUNT m_fcStartOcc;
	FCOUNT m_fcEndOcc;
	SCORE m_scoreGapOpen;
	SCORE m_scoreGapClose;
#if	DOUBLE_AFFINE
	SCORE m_scoreGapOpen2;
	SCORE m_scoreGapClose2;
#endif
//	SCORE m_scoreGapExtend;
	};

struct ProgNode
	{
	ProgNode()
		{
		m_Prof = 0;
		m_EstringL = 0;
		m_EstringR = 0;
		}
	MSA m_MSA;
	ProfPos *m_Prof;
	PWPath m_Path;
	int *m_EstringL;
	int *m_EstringR;
	unsigned m_uLength;
	WEIGHT m_Weight;
	};

extern unsigned ResidueGroup[];
const unsigned RESIDUE_GROUP_MULTIPLE = (unsigned) ~0;

extern PTR_SCOREMATRIX g_ptrScoreMatrix;

ProfPos *ProfileFromMSA(const MSA &a);

SCORE TraceBack(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  PWPath &Path);
SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
void ProgressiveAlign(const SeqVect &v, const Tree &tree, MSA &a);
SCORE MSAPairSP(const MSA &msa1, const MSA &msa2);

void AlignTwoMSAsGivenPath(const PWPath &Path, const MSA &msaA, const MSA &msaB,
  MSA &msaCombined);

void ListProfile(const ProfPos *Prof, unsigned uLength, const MSA *ptrMSA = 0);
SCORE ScoreProfPos2(const ProfPos &PPA, const ProfPos &PPB);
SCORE FastScorePath2(const ProfPos *PA, unsigned uLengthA,
  const ProfPos *PB, unsigned uLengthB, const PWPath &Path);
bool IsHydrophilic(const FCOUNT fcCounts[]);
int PAM200_Letter(unsigned uLetter1, unsigned uLetter2);
SCORE AverageMatchScore(const PWPath &Path, unsigned uEdgeIndex,
  unsigned uWindowLength);
void WindowSmooth(const SCORE Score[], unsigned uCount, unsigned uWindowLength,
  SCORE SmoothScore[], double dCeil = 9e29);
SCORE FastScoreMSA_LA(const MSA &msa, SCORE MatchScore[] = 0);
SCORE FastScoreMSA_NS(const MSA &msa, SCORE MatchScore[] = 0);
SCORE FastScoreMSA_SP(const MSA &msa, SCORE MatchScore[] = 0);
bool RefineMSA(MSA &msa, const Tree &tree);
SCORE MSAQScore(const MSA &msa, SCORE MatchScore[] = 0);
bool RefineBiParts(MSA &msa, const Tree &tree, bool R);
void FindAnchorCols(const MSA &msa, unsigned AnchorCols[],
  unsigned *ptruAnchorColCount);
double PctIdToHeight(double dPctId);
double PctIdToHeightKimura(double dPctId);
double PctIdToHeightMAFFT(double dPctId);
double PctIdToMAFFTDist(double dPctId);
bool RefineBlocks(MSA &msa, const Tree &tree);
bool RefineSubfams(MSA &msaIn, const Tree &tree, unsigned uIters);
void SetMuscleTree(const Tree &tree);
void CalcClustalWWeights(const Tree &tree, WEIGHT Weights[]);
void RealignDiffs(const MSA &msaIn, const Tree &Diffs,
  const unsigned IdToDiffsTreeNodeIndex[], MSA &msaOut);
//void RealignDiffsE(const MSA &msaIn, const SeqVect &v,
//  const Tree &NewTree, const Tree &OldTree,
//  const unsigned uNewNodeIndexToOldNodeIndex[],
//  MSA &msaOut, ProgNode *OldProgNodes);
void RefineTree(MSA &msa, Tree &tree);
bool IsHydrophobic(const FCOUNT fcCounts[]);
void Hydro(ProfPos *Prof, unsigned uLength);
void SetTermGaps(const ProfPos *Prof, unsigned uLength);

// Macros to simulate 2D matrices
#define DPL(PLA, PLB)	DPL_[(PLB)*uPrefixCountA + (PLA)]
#define DPM(PLA, PLB)	DPM_[(PLB)*uPrefixCountA + (PLA)]
#define DPD(PLA, PLB)	DPD_[(PLB)*uPrefixCountA + (PLA)]
#define DPE(PLA, PLB)	DPE_[(PLB)*uPrefixCountA + (PLA)]
#define DPI(PLA, PLB)	DPI_[(PLB)*uPrefixCountA + (PLA)]
#define DPJ(PLA, PLB)	DPJ_[(PLB)*uPrefixCountA + (PLA)]
#define DPU(PLA, PLB)	DPU_[(PLB)*uPrefixCountA + (PLA)]
#define TBM(PLA, PLB)	TBM_[(PLB)*uPrefixCountA + (PLA)]
#define TBD(PLA, PLB)	TBD_[(PLB)*uPrefixCountA + (PLA)]
#define TBE(PLA, PLB)	TBE_[(PLB)*uPrefixCountA + (PLA)]
#define TBI(PLA, PLB)	TBI_[(PLB)*uPrefixCountA + (PLA)]
#define TBJ(PLA, PLB)	TBJ_[(PLB)*uPrefixCountA + (PLA)]

SCORE ScoreProfPos2LA(const ProfPos &PPA, const ProfPos &PPB);
SCORE ScoreProfPos2NS(const ProfPos &PPA, const ProfPos &PPB);
SCORE ScoreProfPos2SP(const ProfPos &PPA, const ProfPos &PPB);
SCORE ScoreProfPos2SPN(const ProfPos &PPA, const ProfPos &PPB);

#endif // FastProf_h
