#include "muscle.h"
#include "profile3.h"
#include "cachemem3.h"
#include "m3alnparams.h"

void BitTraceBack(char **TraceBack, uint uLengthA, uint uLengthB,
  char LastEdge, string &Path);

//static float ScoreProfPos2_LE(const ProfPos3 &PPA, const ProfPos3 &PPB)
//	{
//	float Score = 0;
//	for (unsigned n = 0; n < 20; ++n)
//		{
//		const byte Letter = PPA.m_SortOrder[n];
//		const float FreqA = PPA.m_Freqs[Letter];
//		if (FreqA == 0)
//			break;
//		Score += FreqA*PPB.m_AAScores[Letter];
//		}
//	if (0 == Score)
//		return -2.5;
//	float logScore = logf(Score);
//	float FinalScore = ((logScore - g_CENTER)*(PPA.m_fOcc * PPB.m_fOcc));
//	return FinalScore;
//	}

/***
Should be equivalant:

(a)	-center 2.2 -addcenter 0.0 -normaafreqs 1 -muloccs 1
(b) -center 0.0 -addcenter 2.2 -normaafreqs 0 -muloccs 0

(b) should be faster, saves adding g_CENTER below
***/
float ScoreProfPos2(const ProfPos3 &PPA, const ProfPos3 &PPB)
	{
	//if (g_PPScoreLE)
	//	return ScoreProfPos2_LE(PPA, PPB);

	float Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const byte Letter = PPA.m_SortOrder[n];
		const float FreqA = PPA.m_Freqs[Letter];
		if (FreqA == 0)
			break;
		Score += FreqA*PPB.m_AAScores[Letter];
		}
	//if (g_MulOccs)
	//	Score = (Score + g_CENTER)*(PPA.m_fOcc * PPB.m_fOcc);
	//else
	//	Score = (Score + g_CENTER);
	
	//Score = (Score + g_CENTER);
	return Score;
	}

//static const float MINUS_INFINITY = -9e9f;

#define ALLOC_TRACE()
#define SetDPM(i, j, x)		/* empty  */
#define SetDPD(i, j, x)		/* empty  */
#define SetDPI(i, j, x)		/* empty  */
#define SetTBM(i, j, x)		/* empty  */
#define SetTBD(i, j, x)		/* empty  */
#define SetTBI(i, j, x)		/* empty  */

#define RECURSE_D(i, j)				\
	{								\
	float DD = DRow[j] + e;			\
	float MD = MPrev[j] + PPAs[i-1]->m_GapOpenScore;\
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
	float MI = MCurr[j-1] + PPBs[j-1]->m_GapOpenScore;\
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
	float DM = DRow[j] + PPAs[i-1]->m_GapCloseScore;	\
	float IM = Iij +     PPBs[j-1]->m_GapCloseScore;	\
	float MM = MCurr[j];							\
	TB[i+1][j+1] &= ~BIT_xM;						\
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

static inline void SetBitTBM(char **TB, uint i, uint j, char c)
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
		asserta(false);
		}
	TB[i][j] &= ~BIT_xM;
	TB[i][j] |= Bit;
	}

static inline void SetBitTBD(char **TB, uint i, uint j, char c)
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
		asserta(false);
		}
	TB[i][j] &= ~BIT_xD;
	TB[i][j] |= Bit;
	}

static inline void SetBitTBI(char **TB, uint i, uint j, char c)
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
		asserta(false);
		}
	TB[i][j] &= ~BIT_xI;
	TB[i][j] |= Bit;
	}

float NWSmall3(CacheMem3 &CM, const Profile3 &ProfA,
  const Profile3 &ProfB, string &Path)
	{
	const uint uLengthA = ProfA.GetColCount();
	const uint uLengthB = ProfB.GetColCount();
	const uint uPrefixCountA = uLengthA + 1;
	const uint uPrefixCountB = uLengthB + 1;
	const vector<ProfPos3 *> &PPAs = ProfA.m_PPs;
	const vector<ProfPos3 *> &PPBs = ProfB.m_PPs;

	//const float e = g_GAP_EXT;
	const float e = 0;

	ALLOC_TRACE()

	CM.AllocCache(uPrefixCountA, uPrefixCountB);

	float *MCurr = CM.CacheMCurr;
	float *MNext = CM.CacheMNext;
	float *MPrev = CM.CacheMPrev;
	float *DRow = CM.CacheDRow;

	char **TB = CM.CacheTB;
	for (uint i = 0; i < uPrefixCountA; ++i)
		memset(TB[i], 0, uPrefixCountB);

	float Iij = MINUS_INFINITY;
	SetDPI(0, 0, Iij);

	Iij = PPAs[0]->m_GapOpenScore;
	SetDPI(0, 1, Iij);

	for (uint j = 2; j <= uLengthB; ++j)
		{
		Iij += e;
		SetDPI(0, j, Iij);
		SetTBI(0, j, 'I');
		}

	for (uint j = 0; j <= uLengthB; ++j)
		{
		DRow[j] = MINUS_INFINITY;
		SetDPD(0, j, DRow[j]);
		SetTBD(0, j, 'D');
		}

	MPrev[0] = 0;
	SetDPM(0, 0, MPrev[0]);
	for (uint j = 1; j <= uLengthB; ++j)
		{
		MPrev[j] = MINUS_INFINITY;
		SetDPM(0, j, MPrev[j]);
		}

	MCurr[0] = MINUS_INFINITY;
	SetDPM(1, 0, MCurr[0]);

	MCurr[1] = ScoreProfPos2(*PPAs[0], *PPBs[0]);
	SetDPM(1, 1, MCurr[1]);
	SetBitTBM(TB, 1, 1, 'M');
	SetTBM(1, 1, 'M');

	for (uint j = 2; j <= uLengthB; ++j)
		{
		MCurr[j] = ScoreProfPos2(*PPAs[0], *PPBs[j-1]) +
		  PPBs[0]->m_GapOpenScore + (j - 2)*e + PPBs[j-2]->m_GapCloseScore;
		SetDPM(1, j, MCurr[j]);
		SetBitTBM(TB, 1, j, 'I');
		SetTBM(1, j, 'I');
		}

// Main DP loop
	for (uint i = 1; i < uLengthA; ++i)
		{
		char *TBRow = TB[i];

		Iij = MINUS_INFINITY;
		SetDPI(i, 0, Iij);

		DRow[0] = PPAs[0]->m_GapOpenScore + (i - 1)*e;
		SetDPD(i, 0, DRow[0]);

		MCurr[0] = MINUS_INFINITY; 
		if (i == 1)
			{
			MCurr[1] = ScoreProfPos2(*PPAs[0], *PPBs[0]);
			SetBitTBM(TB, i, 1, 'M');
			SetTBM(i, 1, 'M');
			}
		else
			{
			MCurr[1] = ScoreProfPos2(*PPAs[i-1], *PPBs[0]) +
			  PPAs[0]->m_GapOpenScore + (i - 2)*e + PPAs[i-2]->m_GapCloseScore;
			SetBitTBM(TB, i, 1, 'D');
			SetTBM(i, 1, 'D');
			}
		SetDPM(i, 0, MCurr[0]);
		SetDPM(i, 1, MCurr[1]);

		for (uint j = 1; j < uLengthB; ++j)
			MNext[j+1] = ScoreProfPos2(*PPAs[i], *PPBs[j]);

		for (uint j = 1; j < uLengthB; ++j)
			{
			RECURSE_D(i, j)
			RECURSE_I(i, j)
			RECURSE_M(i, j)
			}
	// Special case for j=uLengthB
		RECURSE_D_BTerm(i)
		RECURSE_I_BTerm(i)

	// Prev := Curr, Curr := Next, Next := Prev
#define Rotate(a, b, c)	{ float *tmp = a; a = b; b = c; c = tmp; }
		Rotate(MPrev, MCurr, MNext);
#undef Rotate
		}

// Special case for i=uLengthA
	char *TBRow = TB[uLengthA];
	MCurr[0] = MINUS_INFINITY;
	if (uLengthA > 1)
		MCurr[1] = ScoreProfPos2(*PPAs[uLengthA-1], *PPBs[0])
		  + (uLengthA - 2)*e +
		  PPAs[0]->m_GapOpenScore + PPAs[uLengthA-2]->m_GapCloseScore;
	else
		MCurr[1] = ScoreProfPos2(*PPAs[uLengthA-1], *PPBs[0]) +
		  PPAs[0]->m_GapOpenScore + PPAs[0]->m_GapCloseScore;
	SetBitTBM(TB, uLengthA, 1, 'D');
	SetTBM(uLengthA, 1, 'D');
	SetDPM(uLengthA, 0, MCurr[0]);
	SetDPM(uLengthA, 1, MCurr[1]);

	DRow[0] = MINUS_INFINITY;
	SetDPD(uLengthA, 0, DRow[0]);
	for (uint j = 1; j <= uLengthB; ++j)
		RECURSE_D_ATerm(j);

	Iij = MINUS_INFINITY;
	for (uint j = 1; j <= uLengthB; ++j)
		RECURSE_I_ATerm(j)

	float MAB = MCurr[uLengthB];
	float DAB = DRow[uLengthB];
	float IAB = Iij;

	float Score = MAB;
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

	BitTraceBack(TB, uLengthA, uLengthB, cEdgeType, Path);

	return Score;
	}
