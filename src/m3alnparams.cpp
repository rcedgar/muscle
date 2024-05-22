#include "muscle.h"
#include "m3alnparams.h"

//M3AlnParams M3AlnParams::m_GlobalM3AlnParams;

void M3AlnParams::Perturb1(float &Param, float MaxDelta) const
	{
	float Sign = (GetRand()%2 == 0 ? -1.0f : +1.0f);

	const uint SMALL_PRIME = 997;
	float f = (GetRand()%SMALL_PRIME)/float(SMALL_PRIME);
	asserta(f >= 0 && f <= 1);
	float Delta = Sign*MaxDelta*f;
	Param += Delta;
	}

void M3AlnParams::Perturb1(float &Param, float MaxDelta)
	{
	float Sign = (GetRand()%2 == 0 ? -1.0f : +1.0f);

	const uint SMALL_PRIME = 997;
	float f = (GetRand()%SMALL_PRIME)/float(SMALL_PRIME);
	asserta(f >= 0 && f <= 1);
	float Delta = Sign*MaxDelta*f;
	Param += Delta;
	}

void M3AlnParams::Print(FILE *f) const
	{
	Log("\n");
	Log("m_GapOpen=%.6g m_Center=%.6g", m_GapOpen, m_Center);
	Log(" linkage=%s", m_Linkage.c_str());
	Log(" treeiters=%u", m_TreeIters);
	Log(" kmerdist=%s\n", m_KmerDist.c_str());
	Log(" perturb(%u)", m_PerturbSeed);
	if (m_PerturbSeed != 0)
		Log(" substmx=%.3g, gapparams=%.3g, distmx=%.3g",
		  m_PerturbSubstMxDelta, m_PerturbGapParamsDelta,
		  m_PerturbDistMxDelta);
	Log("\n");
	for (uint i = 0; i < 3; ++i)
		{
		Log("  SubstMx[%c]: ", g_LetterToChar[i]);
		for (uint j = 0; j < 8; ++j)
			{
			float Score = m_SubstMx_Letter[i][j];
			Log(" %c=%8.4f", g_LetterToChar[j], Score);
			}
		Log("\n");
		}
	Log("\n");

	for (uint i = 0; i < 3; ++i)
		{
		Log("SubstMx-C[%c]: ", g_LetterToChar[i]);
		for (uint j = 0; j < 8; ++j)
			{
			float Score = m_SubstMx_Letter[i][j];
			Log(" %c=%8.4f", g_LetterToChar[j], Score - m_Center);
			}
		Log("\n");
		}
	}

void M3AlnParams::AddCenter(float x)
	{
	if (x == 0)
		return;
	asserta(!m_CenterAdded);
	for (uint i = 0; i < 20; ++i)
		{
		for (uint j = 0; j < 20; ++j)
			m_SubstMx_Letter[i][j] += x;
		}
	m_CenterAdded = true;
	}

void M3AlnParams::SetBlosum(uint PctId, uint n,
  float GapOpen, float Center, uint PerturbSeed, 
  float PerturbSubstMxDelta, float PerturbGapParamsDelta,
  float PerturbDistMxDelta, bool DoLog)
	{
	SetAlphaLC(false);

	m_PerturbSeed = PerturbSeed;
	m_PerturbSubstMxDelta = PerturbSubstMxDelta;
	m_PerturbGapParamsDelta = PerturbGapParamsDelta;
	m_PerturbDistMxDelta = PerturbDistMxDelta;

	m_CenterAdded = false;
	m_PerturbGapParamsDone = false;
	m_PerturbSubstMxDone = false;

	GetSubstMx_Letter_Blosum(PctId, m_SubstMx_Letter);

	float DefaultGapOpen = FLT_MAX;
	float DefaultCenter = FLT_MAX;
	GetGapParams_Blosum(PctId, n, &DefaultGapOpen, &DefaultCenter);

	if (GapOpen != FLT_MAX)
		m_GapOpen = GapOpen;
	else
		m_GapOpen = DefaultGapOpen;

	if (Center != FLT_MAX)
		m_Center = Center;
	else
		m_Center = DefaultCenter;

	AddCenter(m_Center);

	PerturbMyParams();

	m_Ready = true;
	if (DoLog)
		LogMe();
	}

void M3AlnParams::UpdateMx(const Mx2020 &SubstMx_Letter,
  float GapOpen, float Center, bool DoLog)
	{
	SetAlphaLC(false);
	for (uint i = 0; i < 20; ++i)
		for (uint j = 0; j < 20; ++j)
			m_SubstMx_Letter[i][j] = SubstMx_Letter[i][j];

	m_PerturbSubstMxDone = false;
	m_PerturbGapParamsDone = false;
	m_CenterAdded = false;
	m_GapOpen = GapOpen;
	m_Center = Center;
	AddCenter(m_Center);

	m_Ready = true;
	if (DoLog)
		LogMe();
	}

void M3AlnParams::SetFromCmdLine(bool IsNucleo, bool DoLog)
	{
	asserta(!IsNucleo);

	SetAlphaLC(IsNucleo);

	if (optset_substmx)
		{
		asserta(optset_gapopen);
		asserta(optset_center);

		ReadSubstMx_Letter_FromFile(opt(substmx), m_SubstMx_Letter);
		m_GapOpen = (float) opt(gapopen);
		m_Center = (float) opt(center);
		m_PerturbSubstMxDone = false;
		m_PerturbGapParamsDone = false;
		m_CenterAdded = false;
		}
	else
		{
		m_PerturbSubstMxDone = false;
		m_PerturbGapParamsDone = false;
		m_CenterAdded = false;

		uint BlosumPct = optd(blosumpct, 62);
		uint ParamSet = optd(blosumparamset, 0);
		GetSubstMx_Letter_Blosum(BlosumPct, m_SubstMx_Letter);
		GetGapParams_Blosum(BlosumPct, ParamSet, &m_GapOpen, &m_Center);

		if (optset_gapopen)
			m_GapOpen = (float) opt(gapopen);

		if (optset_center)
			m_Center = (float) opt(center);
		}

	AddCenter(m_Center);

	if (optset_perturb)
		{
		m_PerturbSeed = opt(perturb);
		if (m_PerturbSeed != 0)
			{
			m_PerturbSubstMxDelta = 0.1f;
			m_PerturbGapParamsDelta = 0.1f;
			m_PerturbDistMxDelta = 0.1f;
			}
		PerturbMyParams();
		}

	m_Linkage = "biased";
	if (optset_linkage)
		m_Linkage = opt(linkage);

	m_KmerDist = "66";
	if (optset_kmerdist)
		m_KmerDist = opt(kmerdist);

	m_TreeIters = 1;
	if (optset_treeiters)
		m_TreeIters = opt(treeiters);

	m_Ready = true;

	if (DoLog)
		LogMe();
	}

void M3AlnParams::PerturbSubstMx()
	{
	if (m_PerturbSeed == 0 || m_PerturbSubstMxDelta == 0)
		return;
	asserta(!m_PerturbSubstMxDone);
	for (unsigned i = 0; i < 20; ++i)
		for (unsigned j = 0; j < 20; ++j)
			Perturb1(m_SubstMx_Letter[i][j], m_PerturbSubstMxDelta);
	m_PerturbSubstMxDone = true;
	}

void M3AlnParams::PerturbGapParams()
	{
	if (m_PerturbGapParamsDelta == 0)
		return;
	asserta(!m_PerturbGapParamsDone);
	Perturb1(m_GapOpen, m_PerturbGapParamsDelta);
	Perturb1(m_Center, m_PerturbGapParamsDelta);
	m_PerturbGapParamsDone = true;
	}

void M3AlnParams::PerturbMyParams()
	{
	if (m_PerturbSeed == 0)
		return;
	InitPerturb(m_PerturbSeed);
	PerturbGapParams();
	PerturbSubstMx();
	}

void M3AlnParams::InitPerturb(uint Seed)
	{
	m_MinStdRand.seed(Seed);
	}

uint M3AlnParams::GetRand() const
	{
	uint n = m_MinStdRand();
	return n;
	}

void M3AlnParams::PerturbDistMx(vector<vector<float> > &DistMx) const
	{
	if (m_PerturbSeed == 0 || m_PerturbDistMxDelta == 0)
		return;

	const uint N = SIZE(DistMx);
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(DistMx[i]) == N);
		for (uint j = 0; j < i; ++j)
			{
			float d = DistMx[i][j];
			Perturb1(d, m_PerturbDistMxDelta);
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}
