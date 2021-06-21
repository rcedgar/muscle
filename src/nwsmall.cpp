#include "muscle.h"
#include <math.h>
#include "pwpath.h"
#include "profile.h"
#include <stdio.h>

// NW small memory

#define	TRACE	0

#if	TRACE
extern bool g_bKeepSimpleDP;
extern SCORE *g_DPM;
extern SCORE *g_DPD;
extern SCORE *g_DPI;
extern char *g_TBM;
extern char *g_TBD;
extern char *g_TBI;
#endif

#if	TRACE
#define ALLOC_TRACE()								\
	const SCORE UNINIT = MINUS_INFINITY;			\
	const size_t LM = uPrefixCountA*uPrefixCountB;	\
													\
	SCORE *DPM_ = new SCORE[LM];					\
	SCORE *DPD_ = new SCORE[LM];					\
	SCORE *DPI_ = new SCORE[LM];					\
													\
	char *TBM_ = new char[LM];						\
	char *TBD_ = new char[LM];						\
	char *TBI_ = new char[LM];						\
													\
	memset(TBM_, '?', LM);							\
	memset(TBD_, '?', LM);							\
	memset(TBI_, '?', LM);							\
													\
	for (unsigned i = 0; i <= uLengthA; ++i)		\
		for (unsigned j = 0; j <= uLengthB; ++j)	\
			{										\
			DPM(i, j) = UNINIT;						\
			DPD(i, j) = UNINIT;						\
			DPI(i, j) = UNINIT;						\
			}
#else
#define ALLOC_TRACE()
#endif

#if	TRACE
#define SetDPM(i, j, x)		DPM(i, j) = x
#define SetDPD(i, j, x)		DPD(i, j) = x
#define SetDPI(i, j, x)		DPI(i, j) = x
#define SetTBM(i, j, x)		TBM(i, j) = x
#define SetTBD(i, j, x)		TBD(i, j) = x
#define SetTBI(i, j, x)		TBI(i, j) = x
#else
#define SetDPM(i, j, x)		/* empty  */
#define SetDPD(i, j, x)		/* empty  */
#define SetDPI(i, j, x)		/* empty  */
#define SetTBM(i, j, x)		/* empty  */
#define SetTBD(i, j, x)		/* empty  */
#define SetTBI(i, j, x)		/* empty  */
#endif

#define RECURSE_D(i, j)				\
	{								\
	SCORE DD = DRow[j] + e;			\
	SCORE MD = MPrev[j] + PA[i-1].m_scoreGapOpen;\
	if (DD > MD)					\
		{							\
		DRow[j] = DD;				\
		SetTBD(i, j, 'D');			\
		}							\
	else							\
		{							\
		DRow[j] = MD;				\
		/* SetBitTBD(TB, i, j, 'M'); */	\
		TBRow[j] &= ~BIT_xD;		\
		TBRow[j] |= BIT_MD;			\
		SetTBD(i, j, 'M');			\
		}							\
	SetDPD(i, j, DRow[j]);			\
	}

#define RECURSE_D_ATerm(j)	RECURSE_D(uLengthA, j)
#define RECURSE_D_BTerm(j)	RECURSE_D(i, uLengthB)

#define RECURSE_I(i, j)				\
	{								\
	Iij += e;						\
	SCORE MI = MCurr[j-1] + PB[j-1].m_scoreGapOpen;\
	if (MI >= Iij)					\
		{							\
		Iij = MI;					\
		/* SetBitTBI(TB, i, j, 'M'); */	\
		TBRow[j] &= ~BIT_xI;		\
		TBRow[j] |= BIT_MI;			\
		SetTBI(i, j, 'M');			\
		}							\
	else							\
		SetTBI(i, j, 'I');			\
	SetDPI(i, j, Iij);				\
	}

#define RECURSE_I_ATerm(j)	RECURSE_I(uLengthA, j)
#define RECURSE_I_BTerm(j)	RECURSE_I(i, uLengthB)

#define RECURSE_M(i, j)								\
	{												\
	SCORE DM = DRow[j] + PA[i-1].m_scoreGapClose;	\
	SCORE IM = Iij +     PB[j-1].m_scoreGapClose;	\
	SCORE MM = MCurr[j];							\
	TB[i+1][j+1] &= ~BIT_xM;							\
	if (MM >= DM && MM >= IM)						\
		{											\
		MNext[j+1] += MM;							\
		SetDPM(i+1, j+1, MNext[j+1]);				\
		SetTBM(i+1, j+1, 'M');						\
		/* SetBitTBM(TB, i+1, j+1, 'M');	*/		\
		TB[i+1][j+1] |= BIT_MM;						\
		}											\
	else if (DM >= MM && DM >= IM)					\
		{											\
		MNext[j+1] += DM;							\
		SetDPM(i+1, j+1, MNext[j+1]);				\
		SetTBM(i+1, j+1, 'D');						\
		/* SetBitTBM(TB, i+1, j+1, 'D'); */			\
		TB[i+1][j+1] |= BIT_DM;						\
		}											\
	else											\
		{											\
		assert(IM >= MM && IM >= DM);				\
		MNext[j+1] += IM;							\
		SetDPM(i+1, j+1, MNext[j+1]);				\
		SetTBM(i+1, j+1, 'I');						\
		/* SetBitTBM(TB, i+1, j+1, 'I'); */			\
		TB[i+1][j+1] |= BIT_IM;						\
		}											\
	}

#if	TRACE
static bool LocalEq(BASETYPE b1, BASETYPE b2)
	{
	if (b1 < -100000 && b2 < -100000)
		return true;
	double diff = fabs(b1 - b2);
	if (diff < 0.0001)
		return true;
	double sum = fabs(b1) + fabs(b2);
	return diff/sum < 0.005;
	}

static char Get_M_Char(char Bits)
	{
	switch (Bits & BIT_xM)
		{
	case BIT_MM:
		return 'M';
	case BIT_DM:
		return 'D';
	case BIT_IM:
		return 'I';
		}
	Quit("Huh?");
	return '?';
	}

static char Get_D_Char(char Bits)
	{
	return (Bits & BIT_xD) ? 'M' : 'D';
	}

static char Get_I_Char(char Bits)
	{
	return (Bits & BIT_xI) ? 'M' : 'I';
	}

static bool DPEq(char c, SCORE *g_DP, SCORE *DPD_,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	SCORE *DPM_ = g_DP;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			if (!LocalEq(DPM(i, j), DPD(i, j)))
				{
				Log("***DPDIFF*** DP%c(%d, %d) Simple = %.2g, Fast = %.2g\n",
				  c, i, j, DPM(i, j), DPD(i, j));
				return false;
				}
	return true;
	}

static bool CompareTB(char **TB, char *TBM_, char *TBD_, char *TBI_, 
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	SCORE *DPM_ = g_DPM;
	bool Eq = true;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			{
			char c1 = TBM(i, j);
			char c2 = Get_M_Char(TB[i][j]);
			if (c1 != '?' && c1 != c2 && DPM(i, j) > -100000)
				{
				Log("TBM(%d, %d) Simple = %c, NW = %c\n", i, j, c1, c2);
				Eq = false;
				goto D;
				}
			}

D:
	SCORE *DPD_ = g_DPD;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			{
			char c1 = TBD(i, j);
			char c2 = Get_D_Char(TB[i][j]);
			if (c1 != '?' && c1 != c2 && DPD(i, j) > -100000)
				{
				Log("TBD(%d, %d) Simple = %c, NW = %c\n", i, j, c1, c2);
				Eq = false;
				goto I;
				}
			}
I:
	SCORE *DPI_ = g_DPI;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			{
			char c1 = TBI(i, j);
			char c2 = Get_I_Char(TB[i][j]);
			if (c1 != '?' && c1 != c2 && DPI(i, j) > -100000)
				{
				Log("TBI(%d, %d) Simple = %c, NW = %c\n", i, j, c1, c2);
				Eq = false;
				goto Done;
				}
			}
Done:
	if (Eq)
		Log("TB success\n");
	return Eq;
	}

static const char *LocalScoreToStr(SCORE s)
	{
	static char str[16];
	if (s < -100000)
		return "     *";
	sprintf(str, "%6.1f", s);
	return str;
	}

static void LogDP(const SCORE *DPM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %s", LocalScoreToStr(DPM(uPrefixLengthA, uPrefixLengthB)));
		Log("\n");
		}
	}

static void LogBitTB(char **TB, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	Log("Bit TBM:\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = Get_M_Char(TB[uPrefixLengthA][uPrefixLengthB]);
			Log(" %6c", c);
			}
		Log("\n");
		}

	Log("\n");
	Log("Bit TBD:\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = Get_D_Char(TB[uPrefixLengthA][uPrefixLengthB]);
			Log(" %6c", c);
			}
		Log("\n");
		}

	Log("\n");
	Log("Bit TBI:\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = Get_I_Char(TB[uPrefixLengthA][uPrefixLengthB]);
			Log(" %6c", c);
			}
		Log("\n");
		}
	}

static void ListTB(char *TBM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = TBM(uPrefixLengthA, uPrefixLengthB);
			Log(" %6c", c);
			}
		Log("\n");
		}
	}

static const char *BitsToStr(char Bits)
	{
	static char Str[9];

	sprintf(Str, "%cM %cD %cI",
	  Get_M_Char(Bits),
	  Get_D_Char(Bits),
	  Get_I_Char(Bits));
	}
#endif	// TRACE

static inline void SetBitTBM(char **TB, unsigned i, unsigned j, char c)
	{
	char Bit;
	switch (c)
		{
	case 'M':
		Bit = BIT_MM;
		break;
	case 'D':
		Bit = BIT_DM;
		break;
	case 'I':
		Bit = BIT_IM;
		break;
	default:
		Quit("Huh?!");
		}
	TB[i][j] &= ~BIT_xM;
	TB[i][j] |= Bit;
	}

static inline void SetBitTBD(char **TB, unsigned i, unsigned j, char c)
	{
	char Bit;
	switch (c)
		{
	case 'M':
		Bit = BIT_MD;
		break;
	case 'D':
		Bit = BIT_DD;
		break;
	default:
		Quit("Huh?!");
		}
	TB[i][j] &= ~BIT_xD;
	TB[i][j] |= Bit;
	}

static inline void SetBitTBI(char **TB, unsigned i, unsigned j, char c)
	{
	char Bit;
	switch (c)
		{
	case 'M':
		Bit = BIT_MI;
		break;
	case 'I':
		Bit = BIT_II;
		break;
	default:
		Quit("Huh?!");
		}
	TB[i][j] &= ~BIT_xI;
	TB[i][j] |= Bit;
	}

#if	TRACE
#define LogMatrices()											\
	{															\
	Log("Bit DPM:\n");											\
	LogDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPD:\n");											\
	LogDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPI:\n");											\
	LogDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit TB:\n");											\
	LogBitTB(TB, PA, PB, uPrefixCountA, uPrefixCountB);			\
	bool Same;													\
	Same = DPEq('M', g_DPM, DPM_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPM success\n");									\
	Same = DPEq('D', g_DPD, DPD_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPD success\n");									\
	Same = DPEq('I', g_DPI, DPI_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPI success\n");									\
	CompareTB(TB, g_TBM, g_TBD, g_TBI, uPrefixCountA, uPrefixCountB);\
	}
#else
#define LogMatrices()	/* empty */
#endif

static unsigned uCachePrefixCountB;
static unsigned uCachePrefixCountA;
static SCORE *CacheMCurr;
static SCORE *CacheMNext;
static SCORE *CacheMPrev;
static SCORE *CacheDRow;
static char **CacheTB;

static void AllocCache(unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	if (uPrefixCountA <= uCachePrefixCountA && uPrefixCountB <= uCachePrefixCountB)
		return;

	delete[] CacheMCurr;
	delete[] CacheMNext;
	delete[] CacheMPrev;
	delete[] CacheDRow;
	for (unsigned i = 0; i < uCachePrefixCountA; ++i)
		delete[] CacheTB[i];
	delete[] CacheTB;

	uCachePrefixCountA = uPrefixCountA + 1024;
	uCachePrefixCountB = uPrefixCountB + 1024;

	CacheMCurr = new SCORE[uCachePrefixCountB];
	CacheMNext = new SCORE[uCachePrefixCountB];
	CacheMPrev = new SCORE[uCachePrefixCountB];
	CacheDRow = new SCORE[uCachePrefixCountB];

	CacheTB = new char *[uCachePrefixCountA];
	for (unsigned i = 0; i < uCachePrefixCountA; ++i)
		CacheTB[i] = new char [uCachePrefixCountB];
	}

SCORE NWSmall(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (0 == uLengthB || 0 == uLengthA )
		Quit("Internal error, NWSmall: length=0");

	SetTermGaps(PA, uLengthA);
	SetTermGaps(PB, uLengthB);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;
	const SCORE e = g_scoreGapExtend;

	ALLOC_TRACE()

	AllocCache(uPrefixCountA, uPrefixCountB);

	SCORE *MCurr = CacheMCurr;
	SCORE *MNext = CacheMNext;
	SCORE *MPrev = CacheMPrev;
	SCORE *DRow = CacheDRow;

	char **TB = CacheTB;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		memset(TB[i], 0, uPrefixCountB);

	SCORE Iij = MINUS_INFINITY;
	SetDPI(0, 0, Iij);

	Iij = PB[0].m_scoreGapOpen;
	SetDPI(0, 1, Iij);

	for (unsigned j = 2; j <= uLengthB; ++j)
		{
		Iij += e;
		SetDPI(0, j, Iij);
		SetTBI(0, j, 'I');
		}

	for (unsigned j = 0; j <= uLengthB; ++j)
		{
		DRow[j] = MINUS_INFINITY;
		SetDPD(0, j, DRow[j]);
		SetTBD(0, j, 'D');
		}

	MPrev[0] = 0;
	SetDPM(0, 0, MPrev[0]);
	for (unsigned j = 1; j <= uLengthB; ++j)
		{
		MPrev[j] = MINUS_INFINITY;
		SetDPM(0, j, MPrev[j]);
		}

	MCurr[0] = MINUS_INFINITY;
	SetDPM(1, 0, MCurr[0]);

	MCurr[1] = ScoreProfPos2(PA[0], PB[0]);
	SetDPM(1, 1, MCurr[1]);
	SetBitTBM(TB, 1, 1, 'M');
	SetTBM(1, 1, 'M');

	for (unsigned j = 2; j <= uLengthB; ++j)
		{
		MCurr[j] = ScoreProfPos2(PA[0], PB[j-1]) + PB[0].m_scoreGapOpen +
		  (j - 2)*e + PB[j-2].m_scoreGapClose;
		SetDPM(1, j, MCurr[j]);
		SetBitTBM(TB, 1, j, 'I');
		SetTBM(1, j, 'I');
		}

// Main DP loop
	for (unsigned i = 1; i < uLengthA; ++i)
		{
		char *TBRow = TB[i];

		Iij = MINUS_INFINITY;
		SetDPI(i, 0, Iij);

		DRow[0] = PA[0].m_scoreGapOpen + (i - 1)*e;
		SetDPD(i, 0, DRow[0]);

		MCurr[0] = MINUS_INFINITY; 
		if (i == 1)
			{
			MCurr[1] = ScoreProfPos2(PA[0], PB[0]);
			SetBitTBM(TB, i, 1, 'M');
			SetTBM(i, 1, 'M');
			}
		else
			{
			MCurr[1] = ScoreProfPos2(PA[i-1], PB[0]) + PA[0].m_scoreGapOpen +
			  (i - 2)*e + PA[i-2].m_scoreGapClose;
			SetBitTBM(TB, i, 1, 'D');
			SetTBM(i, 1, 'D');
			}
		SetDPM(i, 0, MCurr[0]);
		SetDPM(i, 1, MCurr[1]);

		for (unsigned j = 1; j < uLengthB; ++j)
			MNext[j+1] = ScoreProfPos2(PA[i], PB[j]);

		for (unsigned j = 1; j < uLengthB; ++j)
			{
			RECURSE_D(i, j)
			RECURSE_I(i, j)
			RECURSE_M(i, j)
			}
	// Special case for j=uLengthB
		RECURSE_D_BTerm(i)
		RECURSE_I_BTerm(i)

	// Prev := Curr, Curr := Next, Next := Prev
		Rotate(MPrev, MCurr, MNext);
		}

// Special case for i=uLengthA
	char *TBRow = TB[uLengthA];
	MCurr[0] = MINUS_INFINITY;
	if (uLengthA > 1)
		MCurr[1] = ScoreProfPos2(PA[uLengthA-1], PB[0]) + (uLengthA - 2)*e +
		  PA[0].m_scoreGapOpen + PA[uLengthA-2].m_scoreGapClose;
	else
		MCurr[1] = ScoreProfPos2(PA[uLengthA-1], PB[0]) + PA[0].m_scoreGapOpen +
		  PA[0].m_scoreGapClose;
	SetBitTBM(TB, uLengthA, 1, 'D');
	SetTBM(uLengthA, 1, 'D');
	SetDPM(uLengthA, 0, MCurr[0]);
	SetDPM(uLengthA, 1, MCurr[1]);

	DRow[0] = MINUS_INFINITY;
	SetDPD(uLengthA, 0, DRow[0]);
	for (unsigned j = 1; j <= uLengthB; ++j)
		RECURSE_D_ATerm(j);

	Iij = MINUS_INFINITY;
	for (unsigned j = 1; j <= uLengthB; ++j)
		RECURSE_I_ATerm(j)

	LogMatrices();

	SCORE MAB = MCurr[uLengthB];
	SCORE DAB = DRow[uLengthB];
	SCORE IAB = Iij;

	SCORE Score = MAB;
	char cEdgeType = 'M';
	if (DAB > Score)
		{
		Score = DAB;
		cEdgeType = 'D';
		}
	if (IAB > Score)
		{
		Score = IAB;
		cEdgeType = 'I';
		}

#if TRACE
	Log("    Fast: MAB=%.4g DAB=%.4g IAB=%.4g best=%c\n",
	  MAB, DAB, IAB, cEdgeType);
#endif

	BitTraceBack(TB, uLengthA, uLengthB, cEdgeType, Path);

#if	DBEUG
	Path.Validate();
#endif

	return 0;
	}
