#include "muscle.h"
#include <math.h>
#include "pwpath.h"
#include "profile.h"
#include <stdio.h>

#if	DOUBLE_AFFINE

// NW double affine small memory, term gaps fully penalized
// (so up to caller to adjust in profile if desired).

#define	TRACE	0

#define MIN(x, y)	((x) < (y) ? (x) : (y))

#if	TRACE
extern bool g_bKeepSimpleDP;
extern SCORE *g_DPM;
extern SCORE *g_DPD;
extern SCORE *g_DPE;
extern SCORE *g_DPI;
extern SCORE *g_DPJ;
extern char *g_TBM;
extern char *g_TBD;
extern char *g_TBE;
extern char *g_TBI;
extern char *g_TBJ;
#endif

#if	TRACE
#define ALLOC_TRACE()								\
	const SCORE UNINIT = MINUS_INFINITY;			\
	const size_t LM = uPrefixCountA*uPrefixCountB;	\
													\
	SCORE *DPM_ = new SCORE[LM];					\
	SCORE *DPD_ = new SCORE[LM];					\
	SCORE *DPE_ = new SCORE[LM];					\
	SCORE *DPI_ = new SCORE[LM];					\
	SCORE *DPJ_ = new SCORE[LM];					\
													\
	char *TBM_ = new char[LM];						\
	char *TBD_ = new char[LM];						\
	char *TBE_ = new char[LM];						\
	char *TBI_ = new char[LM];						\
	char *TBJ_ = new char[LM];						\
													\
	memset(TBM_, '?', LM);							\
	memset(TBD_, '?', LM);							\
	memset(TBE_, '?', LM);							\
	memset(TBI_, '?', LM);							\
	memset(TBJ_, '?', LM);							\
													\
	for (unsigned i = 0; i <= uLengthA; ++i)		\
		for (unsigned j = 0; j <= uLengthB; ++j)	\
			{										\
			DPM(i, j) = UNINIT;						\
			DPD(i, j) = UNINIT;						\
			DPE(i, j) = UNINIT;						\
			DPI(i, j) = UNINIT;						\
			DPJ(i, j) = UNINIT;						\
			}
#else
#define ALLOC_TRACE()
#endif

#if	TRACE
#define SetDPM(i, j, x)		DPM(i, j) = x
#define SetDPD(i, j, x)		DPD(i, j) = x
#define SetDPE(i, j, x)		DPE(i, j) = x
#define SetDPI(i, j, x)		DPI(i, j) = x
#define SetDPJ(i, j, x)		DPJ(i, j) = x
#define SetTBM(i, j, x)		TBM(i, j) = x
#define SetTBD(i, j, x)		TBD(i, j) = x
#define SetTBE(i, j, x)		TBE(i, j) = x
#define SetTBI(i, j, x)		TBI(i, j) = x
#define SetTBJ(i, j, x)		TBJ(i, j) = x
#else
#define SetDPM(i, j, x)		/* empty  */
#define SetDPD(i, j, x)		/* empty  */
#define SetDPE(i, j, x)		/* empty  */
#define SetDPI(i, j, x)		/* empty  */
#define SetDPJ(i, j, x)		/* empty  */
#define SetTBM(i, j, x)		/* empty  */
#define SetTBD(i, j, x)		/* empty  */
#define SetTBE(i, j, x)		/* empty  */
#define SetTBI(i, j, x)		/* empty  */
#define SetTBJ(i, j, x)		/* empty  */
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
		SetBitTBD(TB, i, j, 'M');	\
		SetTBD(i, j, 'M');			\
		}							\
	SetDPD(i, j, DRow[j]);			\
	}

#define RECURSE_E(i, j)				\
	{								\
	SCORE EE = ERow[j] + e2;		\
	SCORE ME = MPrev[j] + PA[i-1].m_scoreGapOpen2;\
	if (EE > ME)					\
		{							\
		ERow[j] = EE;				\
		SetTBE(i, j, 'E');			\
		}							\
	else							\
		{							\
		ERow[j] = ME;				\
		SetBitTBE(TB, i, j, 'M');	\
		SetTBE(i, j, 'M');			\
		}							\
	SetDPE(i, j, ERow[j]);			\
	}

#define RECURSE_D_ATerm(j)	RECURSE_D(uLengthA, j)
#define RECURSE_E_ATerm(j)	RECURSE_E(uLengthA, j)

#define RECURSE_D_BTerm(j)	RECURSE_D(i, uLengthB)
#define RECURSE_E_BTerm(j)	RECURSE_E(i, uLengthB)

#define RECURSE_I(i, j)				\
	{								\
	Iij += e;						\
	SCORE MI = MCurr[j-1] + PB[j-1].m_scoreGapOpen;\
	if (MI >= Iij)					\
		{							\
		Iij = MI;					\
		SetBitTBI(TB, i, j, 'M');	\
		SetTBI(i, j, 'M');			\
		}							\
	else							\
		SetTBI(i, j, 'I');			\
	SetDPI(i, j, Iij);				\
	}

#define RECURSE_J(i, j)				\
	{								\
	Jij += e2;						\
	SCORE MJ = MCurr[j-1] + PB[j-1].m_scoreGapOpen2;\
	if (MJ >= Jij)					\
		{							\
		Jij = MJ;					\
		SetBitTBJ(TB, i, j, 'M');	\
		SetTBJ(i, j, 'M');			\
		}							\
	else							\
		SetTBJ(i, j, 'I');			\
	SetDPJ(i, j, Jij);				\
	}

#define RECURSE_I_ATerm(j)	RECURSE_I(uLengthA, j)
#define RECURSE_J_ATerm(j)	RECURSE_J(uLengthA, j)

#define RECURSE_I_BTerm(j)	RECURSE_I(i, uLengthB)
#define RECURSE_J_BTerm(j)	RECURSE_J(i, uLengthB)

#define RECURSE_M(i, j)									\
	{													\
	SCORE Best = MCurr[j];  /*  MM  */					\
	SetTBM(i+1, j+1, 'M');								\
	SetBitTBM(TB, i+1, j+1, 'M');						\
														\
	SCORE DM = DRow[j] + PA[i-1].m_scoreGapClose;		\
	if (DM > Best)										\
		{												\
		Best = DM;										\
		SetTBM(i+1, j+1, 'D');							\
		SetBitTBM(TB, i+1, j+1, 'D');					\
		}												\
														\
	SCORE EM = ERow[j] + PA[i-1].m_scoreGapClose2;		\
	if (EM > Best)										\
		{												\
		Best = EM;										\
		SetTBM(i+1, j+1, 'E');							\
		SetBitTBM(TB, i+1, j+1, 'E');					\
		}												\
														\
	SCORE IM = Iij + PB[j-1].m_scoreGapClose;			\
	if (IM > Best)										\
		{												\
		Best = IM;										\
		SetTBM(i+1, j+1, 'I');							\
		SetBitTBM(TB, i+1, j+1, 'I');					\
		}												\
														\
	SCORE JM = Jij + PB[j-1].m_scoreGapClose2;			\
	if (JM > Best)										\
		{												\
		Best = JM;										\
		SetTBM(i+1, j+1, 'J');							\
		SetBitTBM(TB, i+1, j+1, 'J');					\
		}												\
	MNext[j+1] += Best;									\
	SetDPM(i+1, j+1, MNext[j+1]);						\
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
	case BIT_EM:
		return 'E';
	case BIT_IM:
		return 'I';
	case BIT_JM:
		return 'J';
		}
	Quit("Huh?");
	return '?';
	}

static char Get_D_Char(char Bits)
	{
	return (Bits & BIT_xD) ? 'M' : 'D';
	}

static char Get_E_Char(char Bits)
	{
	return (Bits & BIT_xE) ? 'M' : 'E';
	}

static char Get_I_Char(char Bits)
	{
	return (Bits & BIT_xI) ? 'M' : 'I';
	}

static char Get_J_Char(char Bits)
	{
	return (Bits & BIT_xJ) ? 'M' : 'J';
	}

static bool DPEq(char c, SCORE *g_DP, SCORE *DPD_,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	if (0 == g_DP)
		{
		Log("***DPDIFF*** DP%c=NULL\n", c);
		return true;
		}

	SCORE *DPM_ = g_DP;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			if (!LocalEq(DPM(i, j), DPD(i, j)))
				{
				Log("***DPDIFF*** DP%c(%d, %d) Simple = %.2g, Small = %.2g\n",
				  c, i, j, DPM(i, j), DPD(i, j));
				return false;
				}
	return true;
	}

static bool CompareTB(char **TB, char *TBM_, char *TBD_, char *TBE_, char *TBI_, char *TBJ_,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	if (!g_bKeepSimpleDP)
		return true;
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
				goto E;
				}
			}
E:
	SCORE *DPE_ = g_DPE;
	if (0 == TBE_)
		goto I;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			{
			char c1 = TBE(i, j);
			char c2 = Get_E_Char(TB[i][j]);
			if (c1 != '?' && c1 != c2 && DPE(i, j) > -100000)
				{
				Log("TBE(%d, %d) Simple = %c, NW = %c\n", i, j, c1, c2);
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
				goto J;
				}
			}
J:
	SCORE *DPJ_ = g_DPJ;
	if (0 == DPJ_)
		goto Done;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		for (unsigned j = 0; j < uPrefixCountB; ++j)
			{
			char c1 = TBJ(i, j);
			char c2 = Get_J_Char(TB[i][j]);
			if (c1 != '?' && c1 != c2 && DPJ(i, j) > -100000)
				{
				Log("TBJ(%d, %d) Simple = %c, NW = %c\n", i, j, c1, c2);
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
	Log("Bit TBE:\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = Get_E_Char(TB[uPrefixLengthA][uPrefixLengthB]);
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

	Log("\n");
	Log("Bit TBJ:\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			{
			char c = Get_J_Char(TB[uPrefixLengthA][uPrefixLengthB]);
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
	static char Str[32];

	sprintf(Str, "%cM %cD %cE %cI %cJ",
	  Get_M_Char(Bits),
	  Get_D_Char(Bits),
	  Get_E_Char(Bits),
	  Get_I_Char(Bits),
	  Get_J_Char(Bits));
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
#if	DOUBLE_AFFINE
	case 'E':
		Bit = BIT_EM;
		break;
	case 'I':
		Bit = BIT_IM;
		break;
	case 'J':
		Bit = BIT_JM;
		break;
#endif
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

#if	DOUBLE_AFFINE
static inline void SetBitTBE(char **TB, unsigned i, unsigned j, char c)
	{
	char Bit;
	switch (c)
		{
	case 'M':
		Bit = BIT_ME;
		break;
	case 'E':
		Bit = BIT_EE;
		break;
	default:
		Quit("Huh?!");
		}
	TB[i][j] &= ~BIT_xE;
	TB[i][j] |= Bit;
	}

static inline void SetBitTBJ(char **TB, unsigned i, unsigned j, char c)
	{
	char Bit;
	switch (c)
		{
	case 'M':
		Bit = BIT_MJ;
		break;
	case 'J':
		Bit = BIT_JJ;
		break;
	default:
		Quit("Huh?!");
		}
	TB[i][j] &= ~BIT_xJ;
	TB[i][j] |= Bit;
	}
#endif

#if	TRACE
#define LogMatrices()											\
	{															\
	Log("Bit DPM:\n");											\
	LogDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPD:\n");											\
	LogDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPE:\n");											\
	LogDP(DPE_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPI:\n");											\
	LogDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit DPJ:\n");											\
	LogDP(DPJ_, PA, PB, uPrefixCountA, uPrefixCountB);			\
	Log("Bit TB:\n");											\
	LogBitTB(TB, PA, PB, uPrefixCountA, uPrefixCountB);			\
	bool Same;													\
	Same = DPEq('M', g_DPM, DPM_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPM success\n");									\
	Same = DPEq('D', g_DPD, DPD_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPD success\n");									\
	Same = DPEq('E', g_DPE, DPE_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPE success\n");									\
	Same = DPEq('I', g_DPI, DPI_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPI success\n");									\
	Same = DPEq('J', g_DPJ, DPJ_, uPrefixCountA, uPrefixCountB);\
	if (Same)													\
		Log("DPJ success\n");									\
	CompareTB(TB, g_TBM, g_TBD, g_TBE, g_TBI, g_TBJ, uPrefixCountA, uPrefixCountB);\
	}
#else
#define LogMatrices()	/* empty */
#endif

SCORE NWDASmall(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	ProfPos *pa0 = (ProfPos *) PA;
	ProfPos *pb0 = (ProfPos *) PB;
	ProfPos *paa = (ProfPos *) (PA + uLengthA - 1);
	ProfPos *pbb = (ProfPos *) (PB + uLengthB - 1);

	pa0->m_scoreGapOpen *= -1;
	pb0->m_scoreGapOpen *= -1;

	paa->m_scoreGapClose *= -1;
	pbb->m_scoreGapClose *= -1;

	pa0->m_scoreGapOpen2 *= -1;
	pb0->m_scoreGapOpen2 *= -1;
	paa->m_scoreGapClose2 *= -1;
	pbb->m_scoreGapClose2 *= -1;

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;
	const SCORE e = g_scoreGapExtend;

	const SCORE e2 = g_scoreGapExtend2;
	const SCORE min_e = MIN(g_scoreGapExtend, g_scoreGapExtend2);

	ALLOC_TRACE()

	SCORE *MCurr = new SCORE[uPrefixCountB];
	SCORE *MNext = new SCORE[uPrefixCountB];
	SCORE *MPrev = new SCORE[uPrefixCountB];
	SCORE *DRow = new SCORE[uPrefixCountB];
	SCORE *ERow = new SCORE[uPrefixCountB];

	char **TB = new char *[uPrefixCountA];
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		{
		TB[i] = new char [uPrefixCountB];
		memset(TB[i], 0, uPrefixCountB);
		}

	SCORE Iij = MINUS_INFINITY;
	SetDPI(0, 0, Iij);

	SCORE Jij = MINUS_INFINITY;
	SetDPJ(0, 0, Jij);

	Iij = PB[0].m_scoreGapOpen;
	SetDPI(0, 1, Iij);

	Jij = PB[0].m_scoreGapOpen2;
	SetDPJ(0, 1, Jij);

	for (unsigned j = 2; j <= uLengthB; ++j)
		{
		Iij += e;
		Jij += e2;

		SetDPI(0, j, Iij);
		SetDPJ(0, j, Jij);

		SetTBI(0, j, 'I');
		SetTBJ(0, j, 'J');
		}

	for (unsigned j = 0; j <= uLengthB; ++j)
		{
		DRow[j] = MINUS_INFINITY;
		ERow[j] = MINUS_INFINITY;

		SetDPD(0, j, DRow[j]);
		SetDPE(0, j, ERow[j]);

		SetTBD(0, j, 'D');
		SetTBE(0, j, 'E');
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
		SCORE M = ScoreProfPos2(PA[0], PB[j-1]) + PB[0].m_scoreGapOpen +
		  (j - 2)*e + PB[j-2].m_scoreGapClose;
		SCORE M2 = ScoreProfPos2(PA[0], PB[j-1]) + PB[0].m_scoreGapOpen2 +
		  (j - 2)*e2 + PB[j-2].m_scoreGapClose2;
		
		if (M >= M2)
			{
			MCurr[j] = M;
			SetBitTBM(TB, 1, j, 'I');
			SetTBM(1, j, 'I');
			}
		else
			{
			MCurr[j] = M2;
			SetBitTBM(TB, 1, j, 'J');
			SetTBM(1, j, 'J');
			}
		SetDPM(1, j, MCurr[j]);
		}

// Main DP loop
	for (unsigned i = 1; i < uLengthA; ++i)
		{
		Iij = MINUS_INFINITY;
		Jij = MINUS_INFINITY;
		SetDPI(i, 0, Iij);
		SetDPJ(i, 0, Jij);

		DRow[0] = PA[0].m_scoreGapOpen + (i - 1)*e;
		ERow[0] = PA[0].m_scoreGapOpen2 + (i - 1)*e2;
		SetDPD(i, 0, DRow[0]);
		SetDPE(i, 0, ERow[0]);

		MCurr[0] = MINUS_INFINITY; 
		if (i == 1)
			{
			MCurr[1] = ScoreProfPos2(PA[0], PB[0]);
			SetBitTBM(TB, i, 1, 'M');
			SetTBM(i, 1, 'M');
			}
		else
			{
			SCORE M = ScoreProfPos2(PA[i-1], PB[0]) + PA[0].m_scoreGapOpen +
			  (i - 2)*e + PA[i-2].m_scoreGapClose;
			SCORE M2 = ScoreProfPos2(PA[i-1], PB[0]) + PA[0].m_scoreGapOpen2 +
			  (i - 2)*e2 + PA[i-2].m_scoreGapClose2;
			if (M >= M2)
				{
				MCurr[1] = M;
				SetBitTBM(TB, i, 1, 'D');
				SetTBM(i, 1, 'D');
				}
			else
				{
				MCurr[1] = M2;
				SetBitTBM(TB, i, 1, 'E');
				SetTBM(i, 1, 'E');
				}
			}
		SetDPM(i, 0, MCurr[0]);
		SetDPM(i, 1, MCurr[1]);

		for (unsigned j = 1; j < uLengthB; ++j)
			MNext[j+1] = ScoreProfPos2(PA[i], PB[j]);

		for (unsigned j = 1; j < uLengthB; ++j)
			{
			RECURSE_D(i, j)
			RECURSE_E(i, j)
			RECURSE_I(i, j)
			RECURSE_J(i, j)
			RECURSE_M(i, j)
			}
	// Special case for j=uLengthB
		RECURSE_D_BTerm(i)
		RECURSE_E_BTerm(i)
		RECURSE_I_BTerm(i)
		RECURSE_J_BTerm(i)

	// Prev := Curr, Curr := Next, Next := Prev
		Rotate(MPrev, MCurr, MNext);
		}

// Special case for i=uLengthA
	MCurr[0] = MINUS_INFINITY;
	SCORE M = ScoreProfPos2(PA[uLengthA-1], PB[0]) + (uLengthA - 2)*e +
	  PA[0].m_scoreGapOpen + PA[uLengthA-2].m_scoreGapClose;
	SCORE M2 = ScoreProfPos2(PA[uLengthA-1], PB[0]) + (uLengthA - 2)*e +
	  PA[0].m_scoreGapOpen + PA[uLengthA-2].m_scoreGapClose;
	if (M >= M2)
		{
		MCurr[1] = M;
		SetBitTBM(TB, uLengthA, 1, 'D');
		SetTBM(uLengthA, 1, 'D');
		}
	else
		{
		MCurr[1] = M2;
		SetBitTBM(TB, uLengthA, 1, 'E');
		SetTBM(uLengthA, 1, 'D');
		}
	SetDPM(uLengthA, 0, MCurr[0]);
	SetDPM(uLengthA, 1, MCurr[1]);

	DRow[0] = MINUS_INFINITY;
	ERow[0] = MINUS_INFINITY;

	SetDPD(uLengthA, 0, DRow[0]);
	SetDPE(uLengthA, 0, ERow[0]);

	for (unsigned j = 1; j <= uLengthB; ++j)
		{
		RECURSE_D_ATerm(j);
		RECURSE_E_ATerm(j);
		}

	Iij = MINUS_INFINITY;
	Jij = MINUS_INFINITY;

	for (unsigned j = 1; j <= uLengthB; ++j)
		{
		RECURSE_I_ATerm(j)
		RECURSE_J_ATerm(j)
		}

	LogMatrices();

	SCORE MAB = MCurr[uLengthB];
	SCORE DAB = DRow[uLengthB] + PA[uLengthA-1].m_scoreGapClose;
	SCORE EAB = ERow[uLengthB] + PA[uLengthA-1].m_scoreGapClose2;
	SCORE IAB = Iij + PB[uLengthB-1].m_scoreGapClose;
	SCORE JAB = Jij + PB[uLengthB-1].m_scoreGapClose2;

	SCORE Score = MAB;
	char cEdgeType = 'M';
	if (DAB > Score)
		{
		Score = DAB;
		cEdgeType = 'D';
		}
	if (EAB > Score)
		{
		Score = EAB;
		cEdgeType = 'E';
		}
	if (IAB > Score)
		{
		Score = IAB;
		cEdgeType = 'I';
		}
	if (JAB > Score)
		{
		Score = JAB;
		cEdgeType = 'J';
		}

#if TRACE
	Log("    Small: MAB=%.4g DAB=%.4g EAB=%.4g IAB=%.4g JAB=%.4g best=%c\n",
	  MAB, DAB, EAB, IAB, JAB, cEdgeType);
#endif

	BitTraceBack(TB, uLengthA, uLengthB, cEdgeType, Path);

#if	DBEUG
	Path.Validate();
#endif

	delete[] MCurr;
	delete[] MNext;
	delete[] MPrev;
	delete[] DRow;
	delete[] ERow;
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		delete[] TB[i];
	delete[] TB;

	return 0;
	}
#endif	// DOUBLE_AFFINE
