#ifndef params_h
#define params_h

extern const char *g_pstrInFileName;
extern const char *g_pstrOutFileName;

extern const char *g_pstrFASTAOutFileName;
extern const char *g_pstrMSFOutFileName;
extern const char *g_pstrClwOutFileName;
extern const char *g_pstrClwStrictOutFileName;
extern const char *g_pstrHTMLOutFileName;
extern const char *g_pstrPHYIOutFileName;
extern const char *g_pstrPHYSOutFileName;

extern const char *g_pstrFileName1;
extern const char *g_pstrFileName2;

extern const char *g_pstrSPFileName;
extern const char *g_pstrMatrixFileName;

extern const char *g_pstrUseTreeFileName;
extern bool g_bUseTreeNoWarn;

extern const char *g_pstrComputeWeightsFileName;
extern const char *g_pstrScoreFileName;

extern SCORE g_scoreGapOpen;
extern SCORE g_scoreCenter;
extern SCORE g_scoreGapExtend;
extern SCORE g_scoreGapAmbig;

#if	DOUBLE_AFFINE
extern SCORE g_scoreGapOpen2;
extern SCORE g_scoreGapExtend2;
#endif

extern unsigned g_uSmoothWindowLength;
extern unsigned g_uAnchorSpacing;
extern unsigned g_uMaxTreeRefineIters;

extern unsigned g_uMinDiagLength;
extern unsigned g_uMaxDiagBreak;
extern unsigned g_uDiagMargin;

extern unsigned g_uRefineWindow;
extern unsigned g_uWindowFrom;
extern unsigned g_uWindowTo;
extern unsigned g_uSaveWindow;
extern unsigned g_uWindowOffset;

extern unsigned g_uMaxSubFamCount;

extern unsigned g_uHydrophobicRunLength;
extern float g_dHydroFactor;

extern float g_dSmoothScoreCeil;
extern float g_dMinBestColScore;
extern float g_dMinSmoothScore;
extern float g_dSUEFF;

extern bool g_bPrecompiledCenter;
extern bool g_bNormalizeCounts;
extern bool g_bDiags1;
extern bool g_bDiags2;
extern bool g_bDiags;
extern bool g_bAnchors;
extern bool g_bCatchExceptions;

extern bool g_bMSF;
extern bool g_bAln;
extern bool g_bClwStrict;
extern bool g_bHTML;
extern bool g_bPHYI;
extern bool g_bPHYS;

extern bool g_bQuiet;
extern bool g_bVerbose;
extern bool g_bRefine;
extern bool g_bRefineW;
extern bool g_bRefineX;
extern bool g_bLow;
extern bool g_bSW;
extern bool g_bClusterOnly;
extern bool g_bProfile;
extern bool g_bProfDB;
extern bool g_bPPScore;
extern bool g_bBrenner;
extern bool g_bDimer;
extern bool g_bVersion;
extern bool g_bStable;
extern bool g_bFASTA;
extern bool g_bTomHydro;
extern bool g_bMakeTree;

extern PPSCORE g_PPScore;
extern OBJSCORE g_ObjScore;

extern DISTANCE g_Distance1;
extern CLUSTER g_Cluster1;
extern ROOT g_Root1;
extern SEQWEIGHT g_SeqWeight1;

extern DISTANCE g_Distance2;
extern CLUSTER g_Cluster2;
extern ROOT g_Root2;
extern SEQWEIGHT g_SeqWeight2;

extern unsigned g_uMaxIters;
extern unsigned long g_ulMaxSecs;
extern unsigned g_uMaxMB;

extern SEQTYPE g_SeqType;
extern TERMGAPS g_TermGaps;

#endif // params_h
