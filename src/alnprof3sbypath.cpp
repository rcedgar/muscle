#include "muscle.h"
#include "msa.h"
#include "pwpath.h"
#include "profile3.h"

static ProfPos3 g_PPStart;
static bool InitPPStart()
	{
	g_PPStart.SetStartDimers();
	return true;
	}
static bool g_PPStartInitDone = InitPPStart();

// MM
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetDimersMM(
  const ProfPos3 *PPA, float wA,const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wA*PPA->m_LL + wB*PPB->m_LL;
	PPAB.m_LG = wA*PPA->m_LG + wB*PPB->m_LG;
	PPAB.m_GL = wA*PPA->m_GL + wB*PPB->m_GL;
	PPAB.m_GG = wA*PPA->m_GG + wB*PPB->m_GG;
	}

// MD
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  Bj		(-)
//  X		-	?L	LG
//  -		-	?G	GG
static void SetDimersMD(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wA*PPA->m_LL;
	PPAB.m_LG = wA*PPA->m_LG + wB*(PPB->m_LL + PPB->m_GL);
	PPAB.m_GL = wA*PPA->m_GL;
	PPAB.m_GG = wA*PPA->m_GG + wB*(PPB->m_LG + PPB->m_GG);
	}

// DD
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  (-)		(-)
//  -		-	??	GG
static void SetDimersDD(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wA*PPA->m_LL;
	PPAB.m_LG = wA*PPA->m_LG;
	PPAB.m_GL = wA*PPA->m_GL;
	PPAB.m_GG = wA*PPA->m_GG + wB;
	}

// MI
//  Ai		(-)		Out
//  X		-	?L	LG
//  -		-	?G	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetDimersMI(
  const ProfPos3 *PPA, float wA,const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wB*PPB->m_LL;
	PPAB.m_LG = wB*PPB->m_LG + wA*(PPA->m_LL + PPA->m_GL);
	PPAB.m_GL = wB*PPB->m_GL;
	PPAB.m_GG = wB*PPB->m_GG + wA*(PPA->m_LG + PPA->m_GG);
	}

// DM
//  Ai–1	Ai		Out
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
//  
//  (-)		Bj		
//  -		X		?L	GL
//  -		-		?G	GG
static void SetDimersDM(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wA*PPA->m_LL;
	PPAB.m_LG = wA*PPA->m_LG;
	PPAB.m_GL = wA*PPA->m_GL + wB*(PPB->m_LL + PPB->m_GL);
	PPAB.m_GG = wA*PPA->m_GG + wB*(PPB->m_LG + PPB->m_GG);
	}

// IM
//  (-)		Ai		Out		
//  -		X	?L	GL
//  -		-	?G	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetDimersIM(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wB*PPB->m_LL;
	PPAB.m_LG = wB*PPB->m_LG;
	PPAB.m_GL = wB*PPB->m_GL + wA*(PPA->m_LL + PPA->m_GL);
	PPAB.m_GG = wB*PPB->m_GG + wA*(PPA->m_LG + PPA->m_GG);
	}

// ID
//  (-)		Ai		Out
//  -		X	?L	GL
//  -		-	?G	GG

//  Bj		(-)
//  X		-	?L	LG
//  -		-	?G	GG
static void SetDimersID(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = 0;
	PPAB.m_LG = wB*PPB->m_GL + wB*PPB->m_LL;
	PPAB.m_GL = wA*PPA->m_GL + wA*PPA->m_LL;
	PPAB.m_GG = wA*(PPA->m_LG + PPA->m_GG) + wB*(PPB->m_LG + PPB->m_GG);
	}

// DI
//  Ai		(-)		Out
//  X		-	?L	LG
//  -		-	?G	GG

//  (-)		Bj
//  -		X	?L	GL
//  -		-	?G	GG
static void SetDimersDI(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = 0;
	PPAB.m_LG = wA*PPA->m_GL + wA*PPA->m_LL;
	PPAB.m_GL = wB*PPB->m_GL + wB*PPB->m_LL;
	PPAB.m_GG = wA*(PPA->m_LG + PPA->m_GG) + wB*(PPB->m_LG + PPB->m_GG);
	}

// II
//  (-)		(-)		Out
//  -		-	??	GG

//  Bj–1	Bj
//  X		X	LL	LL
//  X		-	LG	LG
//  -		X	GL	GL
//  -		-	GG	GG
static void SetDimersII(
  const ProfPos3 *PPA, float wA, const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	PPAB.m_LL = wB*PPB->m_LL;
	PPAB.m_LG = wB*PPB->m_LG;
	PPAB.m_GL = wB*PPB->m_GL;
	PPAB.m_GG = wB*PPB->m_GG + wA;
	}

static void SetFreqs1(const ProfPos3 *PP, float w, ProfPos3 &PPAB)
	{
	for (uint i = 0; i < g_AlphaSize; ++i)
		PPAB.m_Freqs[i] = w*PP->m_Freqs[i];
	//PPAB.NormalizeAAFreqsIfRequired();
	}

static void SetFreqs2(const ProfPos3 *PPA, float wA, 
  const ProfPos3 *PPB, float wB, ProfPos3 &PPAB)
	{
	for (uint i = 0; i < g_AlphaSize; ++i)
		PPAB.m_Freqs[i] = wA*PPA->m_Freqs[i] + wB*PPB->m_Freqs[i];
	//PPAB.NormalizeAAFreqsIfRequired();
	}

void AlignTwoProfsGivenPath(const Profile3 &ProfA, float WeightA,
  const Profile3 &ProfB, float WeightB,
  const Mx2020 &SubstMx_Letter, float GapOpen,
  const string &Path, Profile3 &ProfAB)
	{
	asserta(WeightA > 0 && WeightB > 0);
	float wA = WeightA/(WeightA + WeightB);
	float wB = WeightB/(WeightA + WeightB);
	assert(feq(wA + wB, 1.0));

	const uint EdgeCount = SIZE(Path);

#if DEBUG
	{
	uint SumMD = 0;
	uint SumMI = 0;
	for (uint i = 0; i < EdgeCount; ++i)
		{
		switch (Path[i])
			{
		case 'M': ++SumMD; ++SumMI; break;
		case 'D': ++SumMD; break;
		case 'I': ++SumMI; break;
		default:  asserta(false);
			}
		}
	asserta(SumMD == ProfA.GetColCount());
	asserta(SumMI == ProfB.GetColCount());
	}
#endif

	ProfAB.Clear();
	ProfAB.m_PPs.reserve(EdgeCount);
	char cPrevType = 'M';
	uint PosA = 0;
	uint PosB = 0;
	ProfPos3 *PPA = &g_PPStart;
	ProfPos3 *PPB = &g_PPStart;
	for (uint EdgeIndex = 0; EdgeIndex < EdgeCount; ++EdgeIndex)
		{
		const char cType = Path[EdgeIndex];
		ProfPos3 &PPAB = *new ProfPos3;

		PPAB.m_AllGaps = false;
		switch (cType)
			{
		case 'M':
			{
			PPA = ProfA.m_PPs[PosA];
			PPB = ProfB.m_PPs[PosB];
			SetFreqs2(PPA, wA, PPB, wB, PPAB);
			switch (cPrevType)
				{
			case 'M': SetDimersMM(PPA, wA, PPB, wB, PPAB); break;
			case 'D': SetDimersDM(PPA, wA, PPB, wB, PPAB); break;
			case 'I': SetDimersIM(PPA, wA, PPB, wB, PPAB); break;
			default: asserta(false);
				}

			++PosA;
			++PosB;
			break;
			}

		case 'D':
			{
			ProfPos3 *PPA = ProfA.m_PPs[PosA];
			SetFreqs1(PPA, wA, PPAB);
			switch (cPrevType)
				{
			case 'M': SetDimersMD(PPA, wA, PPB, wB, PPAB); break;
			case 'D': SetDimersDD(PPA, wA, PPB, wB, PPAB); break;
			case 'I': SetDimersID(PPA, wA, PPB, wB, PPAB); break;
			default: asserta(false);
				}

			++PosA;
			break;
			}

		case 'I':
			{
			ProfPos3 *PPB = ProfB.m_PPs[PosB];
			SetFreqs1(PPB, wB, PPAB);
			switch (cPrevType)
				{
			case 'M': SetDimersMI(PPA, wA, PPB, wB, PPAB); break;
			case 'D': SetDimersDI(PPA, wA, PPB, wB, PPAB); break;
			case 'I': SetDimersII(PPA, wA, PPB, wB, PPAB); break;
			default: asserta(false);
				}
			++PosB;
			break;
			}

		default:
			assert(false);
			}
		PPAB.SetOcc();
		ProfAB.m_PPs.push_back(&PPAB);
		cPrevType = cType;
		}
	assert(SIZE(ProfAB.m_PPs)== EdgeCount);

	ProfAB.SetScores(SubstMx_Letter, GapOpen);
	ProfAB.Validate();
	}
