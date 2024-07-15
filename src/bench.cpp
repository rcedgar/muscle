#include "muscle.h"
#include "bench.h"
#include "m3alnparams.h"

void Bench::AllocThreads()
	{
	asserta(m_AP != 0);
	if (!m_M3s.empty())
		{
		asserta(SIZE(m_M3s) == m_ThreadCount);
		asserta(SIZE(m_QSs) == m_ThreadCount);
		for (uint i = 0; i < m_ThreadCount; ++i)
			m_M3s[i]->m_AP = m_AP;
		return;
		}
	asserta(m_AP->m_Ready);
	for (uint i = 0; i < m_ThreadCount; ++i)
		{
		Muscle3 *M3 = new Muscle3;
		QScorer *QS = 0;
		QScorer2 *QS2 = 0;
		if (opt(q2))
			QS2 = new QScorer2;
		else
			QS = new QScorer;

		M3->m_AP = m_AP;
		m_M3s.push_back(M3);
		m_QSs.push_back(QS);
		m_QS2s.push_back(QS2);
		}
	}

void Bench::LoadQ2(const string &FileName, const string &FaDir,
  const string &RefDir)
	{
	ReadStringsFromFile(FileName, m_RefNames);

	const uint RefCount = SIZE(m_RefNames);
	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		const string &RefName = m_RefNames[RefIndex];
		const string &RefFileName = RefDir + RefName;

		ProgressStep(RefIndex, RefCount, "Reading %u refs", RefCount);

		MultiSequence *Ref = new MultiSequence;
		Ref->m_DupeLabelsOk = true;
		Ref->LoadMFA(RefFileName, false);
		m_Refs.push_back(Ref);
		}

	for (uint RefIndex = 0; RefIndex < RefCount; ++RefIndex)
		{
		const string &RefName = m_RefNames[RefIndex];
		MultiSequence *MS = new MultiSequence;
		const string &InputFileName = FaDir + RefName;
		MS->LoadMFA(InputFileName, true);
		m_Inputs.push_back(MS);
		}
	}

void Bench::Load(const string &FileName, const string &RefDir)
	{
	if (opt(q2))
		{
		string InDir = opt(indir);
		Dirize(InDir);
		LoadQ2(FileName, InDir, RefDir);
		return;
		}

	ReadStringsFromFile(FileName, m_RefNames);

	const uint MSACount = SIZE(m_RefNames);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &RefName = m_RefNames[MSAIndex];
		const string &RefFileName = RefDir + RefName;

		ProgressStep(MSAIndex, MSACount, "Reading %u MSAs", MSACount);

		MultiSequence *RefMSA = new MultiSequence;
		RefMSA->LoadMFA(RefFileName, false);
		m_Refs.push_back(RefMSA);

		MultiSequence *MS = new MultiSequence;
		MS->LoadMFA(RefFileName, true);
		m_Inputs.push_back(MS);
		}
	}

void Bench::Copy(const Bench &B)
	{
	asserta(m_RefNames.empty());
	asserta(m_Refs.empty());
	asserta(m_Inputs.empty());
	const uint N = SIZE(B.m_RefNames);
	for (uint i = 0; i < N; ++i)
		{
		m_RefNames.push_back(B.m_RefNames[i]);
		m_Refs.push_back(B.m_Refs[i]);
		m_Inputs.push_back(B.m_Inputs[i]);
		}
	}

void Bench::FromSample(const Bench &B, unsigned Pct)
	{
	asserta(m_RefNames.empty());
	asserta(m_Refs.empty());
	asserta(m_Inputs.empty());
	const uint MSACount = SIZE(B.m_RefNames);
	asserta(MSACount > 0);
	asserta(SIZE(B.m_Refs) == MSACount);
	asserta(SIZE(B.m_Inputs) == MSACount);
	uint N = (MSACount*Pct)/100;
	if (N == 0)
		N = 1;
	asserta(N <= MSACount);
	vector<uint> Order;
	for (uint i = 0; i < N; ++i)
		Order.push_back(i);
	Shuffle(Order);
	for (uint i = 0; i < N; ++i)
		{
		uint k = Order[i];
		m_RefNames.push_back(B.m_RefNames[k]);
		m_Refs.push_back(B.m_Refs[k]);
		m_Inputs.push_back(B.m_Inputs[k]);
		}
	}

double Bench::Run(const M3AlnParams &AP)
	{
	m_AP = &AP;
	m_FinalScore = DBL_MAX;
	m_ThreadCount = GetRequestedThreadCount();
	AllocThreads();
	m_TCs.clear();

	const uint MSACount = SIZE(m_Inputs);
	asserta(MSACount > 0);
	asserta(SIZE(m_Refs) == MSACount);
	m_TCs.resize(MSACount, DBL_MAX);
	uint Counter = 0;
	double SumQ = 0;
	double SumTC = 0;
	time_t LastProgressTime = time(0);
	time_t t = 0;
	uint SharedMSAIndex = 0;
#pragma omp parallel for num_threads(m_ThreadCount)
	for (int k = 0; k < (int) MSACount; ++k)
		{
		uint MSAIndex = UINT_MAX;
#pragma omp critical
		{
		MSAIndex = SharedMSAIndex;
		++SharedMSAIndex;
		}

		asserta(MSAIndex < MSACount);
		uint ThreadIndex = GetThreadIndex();
		const string &RefName = m_RefNames[MSAIndex];
		Muscle3 *M3 = m_M3s[ThreadIndex];
		QScorer *QS = m_QSs[ThreadIndex];
		QScorer2 *QS2 = m_QS2s[ThreadIndex];
		const MultiSequence *Input = m_Inputs[MSAIndex];
		const MultiSequence *Ref = m_Refs[MSAIndex];
		M3->Run(AP, *Input);
		const MultiSequence *TestMSA = M3->m_FinalMSA;
		double Q = 0;
		double TC = 0;
		if (opt(q2))
			{
			Q = QS2->Run(*TestMSA, *Ref);
			TC = Q;
			}
		else
			{
			QS->Run(RefName, *TestMSA, *Ref);
			Q = QS->m_Q;
			TC = QS->m_TC;
			}
#pragma omp critical
		{
		m_TCs[MSAIndex] = TC;

		++Counter;
		SumQ += Q;
		SumTC += TC;
		m_MeanQ = SumQ/Counter;
		m_MeanTC = SumTC/Counter;

		if (m_ShowProgress && ((t = time(0)) - LastProgressTime >= 1))
			{
			const char *Prefix = GetProgressPrefixCStr();
			printf("%s  AvgQ=%.3f AvgTC=%.3f [%u / %u]       \r",
			  Prefix, m_MeanQ, m_MeanTC, Counter, MSACount);
			LastProgressTime = t;
			}
		}
		}
	if (m_ShowProgress)
		printf("\n");

	m_MeanQ = SumQ/MSACount;
	m_MeanTC = SumTC/MSACount;
	if (m_ShowProgress)
		printf("  AvgQ=%6.4f AvgTC=%6.4f N=%u\n",
			m_MeanQ, m_MeanTC, MSACount);
	m_FinalScore = m_MeanQ + m_MeanTC;
	return m_FinalScore;
	}

void Bench::TCsToFile(const string &FileName) const
	{
	if (FileName.empty())
		return;
	const uint MSACount = SIZE(m_RefNames);
	asserta(SIZE(m_TCs) == MSACount);
	FILE *f = CreateStdioFile(FileName);
	for (uint i = 0; i < MSACount; ++i)
		{
		const char *RefName = m_RefNames[i].c_str();
		fprintf(f, "%s\t%.4f\n", RefName, m_TCs[i]);
		}
	CloseStdioFile(f);
	}
