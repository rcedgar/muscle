#include "muscle.h"
#include "sweeper.h"
#include "sort.h"

float Sweeper::GetRandomValue(float Lo, float Hi) const
	{
	asserta(Hi > Lo);
	const uint SMALL_PRIME = 997;
	uint r = randu32()%SMALL_PRIME;
	float f = r/float(SMALL_PRIME);
	asserta(f >= 0.0f && f <= 1.0f);
	float Value = Lo + f*(Hi - Lo);
	asserta(Value >= Lo && Value <= Hi);
	return Value;
	}

void Sweeper::GetRandomValues(const vector<float> &Los,
  const vector<float> &His, vector<float> &Values) const
	{
	Values.clear();
	for (uint i = 0; i < m_ParamCount; ++i)
		{
		float Value = GetRandomValue(Los[i], His[i]);
		Values.push_back(Value);
		}
	}

float Sweeper::GetGridValue(float Lo, float Hi, uint Step, uint n) const
	{
	asserta(Lo < Hi);
	asserta(n > 0);
	float dStep = float(Step);
	if (m_GridNoiseFract != 0)
		{
		asserta(m_GridNoiseFract > 0 && m_GridNoiseFract < 1);
		float d = GetRandomValue(-m_GridNoiseFract, m_GridNoiseFract);
		dStep += d;
		}

	float Value = FLT_MAX;
	if (n == 1)
		Value = Lo;
	else
		Value = Lo + (Hi - Lo)*dStep/(n-1);
	return Value;
	}

void Sweeper::ExploreGrid(const vector<float> &Los,
  const vector<float> &His, const vector<uint> &Sizes)
	{
	asserta(SIZE(Los) == m_ParamCount);
	asserta(SIZE(His) == m_ParamCount);
	asserta(SIZE(Sizes) == m_ParamCount);

	m_GridLos = &Los;
	m_GridHis = &His;
	m_GridSizes = &Sizes;

	m_GridCount = Sizes[0];
	for (uint i = 1; i < m_ParamCount; ++i)
		m_GridCount *= Sizes[i];

	m_GridCoords.clear();
	m_GridCoords.resize(m_ParamCount, 0);
	for (m_GridCounter = 0; m_GridCounter < m_GridCount; ++m_GridCounter)
		{
		vector<float> Values;
		for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
			{
			float Lo = Los[ParamIndex];
			float Hi = His[ParamIndex];
			uint n = Sizes[ParamIndex];
			uint Step = m_GridCoords[ParamIndex];
			float Value = GetGridValue(Lo, Hi, Step, n);
			Values.push_back(Value);
			}
		Run1(Values);
		bool Ok = GetNextEnumGrid(Sizes, m_GridCoords);
		if (m_GridCounter + 1 == m_GridCount)
			asserta(!Ok);
		else
			asserta(Ok);
		}
	m_GridCounter = UINT_MAX;
	m_GridCount = UINT_MAX;
	}

void Sweeper::ExploreRandom(const vector<float> &Los,
  const vector<float> &His, uint n)
	{
	for (uint i = 0; i < n; ++i)
		{
		vector<float> Values;
		GetRandomValues(Los, His, Values);
		Run1(Values);
		}
	}

void Sweeper::SetFev(const string &FileName)
	{
	asserta(m_fFev == 0);
	m_fFev = CreateStdioFile(FileName);
	}

void Sweeper::SetParamNames(const vector<string> &ParamNames)
	{
	m_ParamNames = ParamNames;
	m_ParamCount = SIZE(m_ParamNames);
	}

void Sweeper::Run1(const vector<float> &ParamValues)
	{
	m_ParamValuesVec.push_back(ParamValues);
	double Q, TC;
	(m_GetScore)(*this, ParamValues, Q, TC);
	//float Score = float(Q + TC)/2.0f;
	float Score = (float) TC;
	m_Scores.push_back(Score);
	m_Qs.push_back(float(Q));
	m_TCs.push_back(float(TC));
	uint Index = SIZE(m_Scores);
	asserta(SIZE(m_ParamValuesVec) == Index);
	bool NewBest = false;
	if (Score > m_BestScore)
		{
		m_BestScore = Score;
		m_BestIndexes.clear();
		m_BestIndexes.push_back(Index);
		NewBest = true;
		}
	else if (Score == m_BestScore)
		{
		m_BestIndexes.push_back(Index);
		NewBest = true;
		}

	if (m_fFev != 0)
		{
		fprintf(m_fFev, "%u", Index);
		fprintf(m_fFev, "\tscore=%.8g", Score);
		fprintf(m_fFev, "\tQ=%.8g", Q);
		fprintf(m_fFev, "\tTC=%.8g", TC);
		for (uint i = 0; i < m_ParamCount; ++i)
			fprintf(m_fFev, "\t%s=%.8g",
			  m_ParamNames[i].c_str(), ParamValues[i]);
		if (!m_GridCoords.empty())
			{
			asserta(SIZE(m_GridCoords) == m_ParamCount);
			asserta(m_GridSizes != 0);
			asserta(SIZE(*m_GridSizes) == m_ParamCount);
			for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
				{
				uint Coord = m_GridCoords[ParamIndex];
				uint Size = (*m_GridSizes)[ParamIndex];
				fprintf(m_fFev, "\tgridcoord%u=%u/%u",
				  ParamIndex, Coord, Size);
				}
			}
		if (NewBest)
			fprintf(m_fFev, "\tnewbest=yes");
		fprintf(m_fFev, "\n");
		fflush(m_fFev);
		}
	}

void Sweeper::GetSortOrder(vector<uint> &Order) const
	{
	const uint N = SIZE(m_Scores);
	asserta(SIZE(m_ParamValuesVec) == N);

	uint *Indexes = myalloc(uint, N);
	for (uint i = 0; i < N; ++i)
		Indexes[i] = i;

	Order.resize(N);
	QuickSortOrderDesc(m_Scores.data(), N, Order.data());
	}

void Sweeper::LogDistinctTop(const vector<float> &Deltas, uint MaxCloudSize,
  float MaxCloudScoreDist, float MaxCloudParamDist, uint n) const
	{
	vector<uint> Indexes;
	GetDistinctTopIndexes(n, m_MaxDistinctScoreDrop, m_MinDistinctParamDist,
	  Indexes);
	const uint NIX = SIZE(Indexes);
	uint M = min(NIX, n);
	Log("Distinct top (%u of %u)\n", M, NIX);

	Log(" n");
	Log("  %5.5s", "Cloud");
	Log("  %8.8s", "Q");
	Log("  %8.8s", "TC");
	Log("  %8.8s", "ParamDst");
	Log("  %8.8s", "ScoreDff");
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		Log("  %12.12s", m_ParamNames[ParamIndex].c_str());
	Log("\n");

	for (uint Rank = 0; Rank < M; ++Rank)
		{
		uint Index = Indexes[Rank];
		const vector<float> &Params = GetParams(Index);

		float Score = GetScore(Index);
		float Q = GetQ(Index);
		float TC = GetTC(Index);

		Log("%2u", Rank);
		Log("  %5.5s", "");
		Log("  %8.5f", GetQ(Index));
		Log("  %8.5f", GetTC(Index));
		Log("  %8.8s", "");
		Log("  %8.8s", "");
		for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
			Log("  %12.8g", Params[ParamIndex]);
		Log("\n");

		vector<uint> CloudIndexes;
		GetCloud(Index, MaxCloudSize, MaxCloudScoreDist, MaxCloudParamDist,
		  CloudIndexes);
		uint CN = SIZE(CloudIndexes);
		for (uint ci = 0; ci < CN; ++ci)
			{
			uint Index2 = CloudIndexes[ci];
			const vector<float> &Params2 = GetParams(Index2);
			float Score2 = GetScore(Index2);
			float Q2 = GetQ(Index2);
			float TC2 = GetTC(Index2);
			float ParamDist = GetParamDist(Index, Index2);
			float ScoreDiff = GetScoreDiff(Index, Index2);
			
			Log("%2u", Rank);
			Log("  %5u", ci);
			Log("  %8.5f", Q2);
			Log("  %8.5f", TC2);
			Log("  %8.5f", ParamDist);
			Log("  %8.5f", ScoreDiff);
			for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
				Log("  %12.8g", Params2[ParamIndex]);
			Log("\n");
			}
		if (CN > 0)
			Log("\n");
		}
	}

void Sweeper::LogIndexes(const vector<uint> &Indexes) const
	{
	const uint N = SIZE(Indexes);
	Log(" Score");
	Log("       Q");
	Log("      TC");
	for (uint i = 0; i < m_ParamCount; ++i)
		Log("  %12.12s", m_ParamNames[i].c_str());
	Log("\n");
	for (uint k = 0; k < SIZE(Indexes); ++k)
		{
		uint i = Indexes[k];
		float Score = m_Scores[i];
		float Q = m_Qs[i];
		float TC = m_TCs[i];
		const vector<float> &Values = m_ParamValuesVec[i];

		Log("%6.4f", Score);
		Log("  %6.4f", Q);
		Log("  %6.4f", TC);
		for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
			{
			float Value = m_ParamValuesVec[i][ParamIndex];
			Log("%10.4g", Value);
			}
		Log("\n");
		}
	}

void Sweeper::LogTop(uint n) const
	{
	const uint N = SIZE(m_Scores);
	vector<uint> Order;
	GetSortOrder(Order);
	uint M = min(N, n);
	Log("Top params:\n");
	Log(" Score");
	Log("       Q");
	Log("      TC");
	for (uint i = 0; i < m_ParamCount; ++i)
		Log("  %8.8s", m_ParamNames[i].c_str());
	Log("\n");
	for (uint k = 0; k < M; ++k)
		{
		uint i = Order[k];
		float Score = m_Scores[i];
		float Q = m_Qs[i];
		float TC = m_TCs[i];
		const vector<float> &Values = m_ParamValuesVec[i];

		Log("%6.4f", Score);
		Log("  %6.4f", Q);
		Log("  %6.4f", TC);
		for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
			{
			float Value = m_ParamValuesVec[i][ParamIndex];
			Log("%12.8g", Value);
			}
		Log("\n");
		}
	}

void Sweeper::GetSpatteredValues(const vector<float> &CenterValues,
  const vector<float> &MaxDeltas, vector<float> &Values) const
	{
	Values.clear();
	asserta(SIZE(CenterValues) == m_ParamCount);
	asserta(SIZE(MaxDeltas) == m_ParamCount);
	Values.clear();
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		{
		float CenterValue = CenterValues[ParamIndex];
		float MaxDelta = MaxDeltas[ParamIndex];
		const uint M = 1000003;
		uint r = randu32()%M;
		float f = float(r)/M;
		asserta(f >= 0.0f && f <= 1.0f);
		float g = 0.2f + f*0.8f;
		float Delta = MaxDelta*g;
		if (r%2 == 0)
			Delta = -Delta;
		float Value = CenterValue + Delta;
		Values.push_back(Value);
		}
	}

// Maximize parameter diversity
void Sweeper::GetDistinctTopIndexes(uint MaxCount, float MaxScoreDrop,
  float MinParamDist, vector<uint> &Indexes) const
	{
	Indexes.clear();
	const uint N = SIZE(m_Scores);
	if (N == 0)
		return;
	vector<uint> Order;
	GetSortOrder(Order);
	asserta(SIZE(Order) == N);
	uint TopIndex = Order[0];
	Indexes.push_back(TopIndex);
	for (uint k = 0; k < N; ++k)
		{
		uint Index = Order[k];
		float ScoreDiff = GetScoreDiff(TopIndex, Index);
		asserta(ScoreDiff >= 0);
		if (ScoreDiff > MaxScoreDrop)
			{
			Order.resize(k+1);
			break;
			}
		}

// Now Order[] has all indexes with high enough scores
	uint M = SIZE(Order);
	vector<float> ParamDists;
	for (uint i = 0; i < M; ++i)
		{
		uint Index = Order[i];
		bool RejectParamsTooSimilar = false;
		for (uint j = 0; j < SIZE(Indexes); ++j)
			{
			float Dist = GetParamDist(Index, Indexes[j]);
			if (Dist < MinParamDist)
				{
				RejectParamsTooSimilar = true;
				break;
				}
			}
		if (RejectParamsTooSimilar)
			continue;
		Indexes.push_back(Index);
		if (SIZE(Indexes) >= MaxCount)
			return;
		}
	}

void Sweeper::ExploreSpatter(const vector<vector<float> > &StartValueVec,
  const vector<float> &StartDeltas, uint TriesPerIter,
  uint MaxIters, uint MaxFailIters, float Shrink)
	{
	const uint StartCount = SIZE(StartValueVec);
	asserta(StartCount > 0);
	for (uint StartIndex = 0; StartIndex < StartCount; ++StartIndex)
		{
		const vector<float> &StartValues = StartValueVec[StartIndex];
		asserta(SIZE(StartValues) == m_ParamCount);
		asserta(SIZE(StartDeltas) == m_ParamCount);
		Run1(StartValues);
		}
	asserta(SIZE(m_Scores) == StartCount);
	asserta(SIZE(m_ParamValuesVec) == StartCount);

	m_SpatterSeedIndexes.clear();
	m_SpatterSeedIndexes.push_back(0);

	m_SpatterTriesPerIter = TriesPerIter;

	asserta(Shrink > 0 && Shrink < 1);
	m_SpatterShrink = Shrink;

	m_SpatterDeltas = StartDeltas;

	m_MaxDistinctScoreDrop = m_StartMaxDistinctScoreDrop;
	m_MinDistinctParamDist = m_StartMinDistinctParamDist;

	for (m_SpatterIter = 0; m_SpatterIter < MaxIters; ++m_SpatterIter)
		{
///////////////////////////////////////////////////////////
		if (m_SpatterIter == 0)
			m_SpatterTriesPerIter = 4*TriesPerIter;
		else
			m_SpatterTriesPerIter = TriesPerIter;
///////////////////////////////////////////////////////////

		printf("\nIter %u failed %u/%u", m_SpatterIter + 1,
		  m_SpatterFailedIterCount, MaxFailIters);
		printf(" scoredrop=%.3g paramdist=%.3g",
		  m_MaxDistinctScoreDrop, m_MinDistinctParamDist);
		PrintfSpatterDeltas();

		Log("\n");
		Log(">>>>> m_SpatterIter=%u (max %u, failed %u / %u) <<<<<\n",
		  m_SpatterIter, MaxIters, m_SpatterFailedIterCount, MaxFailIters);
		Log("m_MaxDistinctScoreDrop = %.4g, m_MinDistinctParamDist = %.4g\n",
		  m_MaxDistinctScoreDrop, m_MinDistinctParamDist);
		LogSpatterDeltas();

///////////////////////////////////////////////////////////
		bool Improved = SpatterIter();
///////////////////////////////////////////////////////////
		m_SpatterTriesPerIter = TriesPerIter;

		Log(">>>>> Improved = %c\n", yon(Improved)),

		LogTop(10);
		LogDistinctTop(m_SpatterDeltas, 8,
		  m_MaxDistinctScoreDrop, m_MinDistinctParamDist, 10);

		if (Improved)
			m_SpatterFailedIterCount = 0;
		else
			{
			++m_SpatterFailedIterCount;
			if (m_SpatterFailedIterCount >= MaxFailIters)
				{
				printf("\nConverged, max failed iters\n");
				Log("\nConverged, max failed iters\n");
				break;
				}
			}
		SpatterUpdateDeltas_Shrink(m_SpatterShrink);

		if (m_MaxDistinctScoreDrop*0.9f >= m_EndMaxDistinctScoreDrop)
			m_MaxDistinctScoreDrop *= 0.9f;
		if (m_MinDistinctParamDist*0.9f >= m_EndMinDistinctParamDist)
			m_MinDistinctParamDist *= 0.8f;
		}
	}

void Sweeper::SpatterUpdateDeltas_Shrink(float Shrink)
	{
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		{
		if (m_SpatterDeltas[ParamIndex] > m_MinDelta)
			m_SpatterDeltas[ParamIndex] *= Shrink;
		}
	}

void Sweeper::PrintfSpatterDeltas() const
	{
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		{
		printf("  %s", m_ParamNames[ParamIndex].c_str());
		printf("+%.4g", m_SpatterDeltas[ParamIndex]);
		}
	printf("\n");
	}

void Sweeper::LogSpatterDeltas() const
	{
	Log("Deltas:");
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		{
		Log("  %s", m_ParamNames[ParamIndex].c_str());
		Log("=%.4g", m_SpatterDeltas[ParamIndex]);
		}
	Log("\n");
	}

bool Sweeper::SpatterIter()
	{
	vector<uint> NewSeedIndexes;
	GetDistinctTopIndexes(m_SpatterTriesPerIter,
	  m_MaxDistinctScoreDrop, m_MinDistinctParamDist,
	  NewSeedIndexes);
	asserta(!NewSeedIndexes.empty());
	m_SpatterSeedIndexes = NewSeedIndexes;

	Log("\nIter %u seeds (%u):\n",
	  m_SpatterIter, SIZE(m_SpatterSeedIndexes));
	LogIndexes(m_SpatterSeedIndexes);

	float BestScoreStart = m_BestScore;
	uint K = SIZE(m_SpatterSeedIndexes);
	asserta(K > 0);
	uint k = 0;
	for (m_SpatterTry = 0; m_SpatterTry < m_SpatterTriesPerIter; ++m_SpatterTry)
		{
		asserta(k < K);
		uint Index = m_SpatterSeedIndexes[k];
		const vector<float> &Values = m_ParamValuesVec[Index];
		k = k%K;

		vector<float> SpatteredValues;
		GetSpatteredValues(Values, m_SpatterDeltas, SpatteredValues);

		Run1(SpatteredValues);
		}

	float BestScoreEnd = m_BestScore;
	bool Improved = (BestScoreEnd - BestScoreStart) > 0.001;
	return Improved;
	}

const vector<float> &Sweeper::GetParams(uint Index) const
	{
	asserta(Index < SIZE(m_ParamValuesVec));
	return m_ParamValuesVec[Index];
	}

float Sweeper::GetParamDist(uint Index1, uint Index2) const
	{
	const vector<float> &Params1 = GetParams(Index1);
	const vector<float> &Params2 = GetParams(Index2);
	float Sum2 = 0;
	for (uint ParamIndex = 0; ParamIndex < m_ParamCount; ++ParamIndex)
		{
		float Param1 = Params1[ParamIndex];
		float Param2 = Params2[ParamIndex];
		float d = Param1 - Param2;
		Sum2 += d*d;
		}
	float Dist = sqrt(Sum2);
	return Dist;
	}

float Sweeper::GetScore(uint Index) const
	{
	asserta(Index < SIZE(m_Scores));
	float Score = m_Scores[Index];
	return Score;
	}

float Sweeper::GetQ(uint Index) const
	{
	asserta(Index < SIZE(m_Qs));
	float Score = m_Qs[Index];
	return Score;
	}

float Sweeper::GetTC(uint Index) const
	{
	asserta(Index < SIZE(m_TCs));
	float TC = m_TCs[Index];
	return TC;
	}

float Sweeper::GetScoreDiff(uint Index1, uint Index2) const
	{
	float Score1 = GetScore(Index1);
	float Score2 = GetScore(Index2);
	return Score1 - Score2;
	}

float Sweeper::GetScoreDist(uint Index1, uint Index2) const
	{
	float Diff = GetScoreDiff(Index1, Index2);
	return fabs(Diff);
	}

void Sweeper::GetCloud(uint Index, uint MaxSize, float MaxScoreDist,
  float MaxParamDist, vector<uint> &Indexes) const
	{
	Indexes.clear();
	const uint N = SIZE(m_ParamValuesVec);
	vector<uint> TmpIndexes;
	vector<float> ScoreDiffs;
	for (uint Index2 = 0; Index2 < N; ++Index2)
		{
		if (Index2 == Index)
			continue;
		const vector<float> &Params = GetParams(Index);
		float ScoreDiff = GetScoreDiff(Index, Index2);
		if (fabs(ScoreDiff) > MaxScoreDist)
			continue;
		float ParamDist = GetParamDist(Index, Index2);
		if (ParamDist > MaxParamDist)
			continue;
		TmpIndexes.push_back(Index2);
		ScoreDiffs.push_back(ScoreDiff);
		}
	const uint M = SIZE(TmpIndexes);
	if (M == 0)
		return;
	vector<uint> Order(M);
	QuickSortOrder(ScoreDiffs.data(), M, Order.data());
	for (uint i = 0; i < min(M, MaxSize); ++i)
		Indexes.push_back(TmpIndexes[Order[i]]);
	}
