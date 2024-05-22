#pragma once

class Sweeper;

typedef void (*ptrfn_SweeperGetScore)
  (const Sweeper &S, const vector<float> &ParamValues,
  double &Q, double &TC);

class Sweeper
	{
public:
	vector<string> m_ParamNames;
	vector<float> m_Scores;
	vector<float> m_Qs;
	vector<float> m_TCs;
	vector<vector<float> > m_ParamValuesVec;
	float m_BestScore = 0;
	vector<uint> m_BestIndexes;
	ptrfn_SweeperGetScore m_GetScore = 0;
	float m_GridNoiseFract = 0;

	uint m_ParamCount = 0;

	FILE *m_fFev = 0;

// Grid search
	const vector<float> *m_GridLos = 0;
	const vector<float> *m_GridHis = 0;
	const vector<uint> *m_GridSizes = 0;
	vector<uint> m_GridCoords;
	uint m_GridCounter = UINT_MAX;
	uint m_GridCount = UINT_MAX;

// Spatter search
	uint m_SpatterTriesPerIter = UINT_MAX;
	uint m_SpatterTry = UINT_MAX;
	uint m_SpatterIter = UINT_MAX;
	float m_SpatterShrink = FLT_MAX;
	float m_MinDelta = 0.05f;
	uint m_SpatterFailedIterCount = 0;
	vector<uint> m_SpatterSeedIndexes;
	vector<float> m_SpatterDeltas;

	float m_StartMaxDistinctScoreDrop = 0.04f;
	float m_EndMaxDistinctScoreDrop = 0.01f;
	float m_StartMinDistinctParamDist = 1.0f;
	float m_EndMinDistinctParamDist = 0.2f;

	float m_MaxDistinctScoreDrop = FLT_MAX;
	float m_MinDistinctParamDist = FLT_MAX;

public:
	void SetParamNames(const vector<string> &ParamNames);
	void SetFev(const string &FevFileName);
	void Run1(const vector<float> &ParamValues);
	void GetRandomValues(const vector<float> &Los,
	  const vector<float> &His, vector<float> &Values) const;
	float GetRandomValue(float Lo, float Hi) const;
	void ExploreRandom(const vector<float> &Los,
	  const vector<float> &His, uint n);
	void GetSpatteredValues(const vector<float> &CenterValues,
	  const vector<float> &MaxDeltas, vector<float> &Values) const;
	float GetGridValue(float Lo, float Hi, uint Step, uint n) const;
	void ExploreGrid(const vector<float> &Los,
	  const vector<float> &His, const vector<uint> &ns);
	void LogTop(uint n = 10) const;
	void LogDistinctTop(const vector<float> &Deltas, uint MaxCloudSize,
	  float MaxCloudScoreDist, float MaxCloudParamDist, uint n) const;
	void GetSortOrder(vector<uint> &Order) const;
	void GetDistinctTopIndexes(uint MaxCount, float MaxScoreDrop,
	  float MinParamDist, vector<uint> &Indexes) const;
	void LogIndexes(const vector<uint> &Indexes) const;
	void ExploreSpatter(const vector<vector<float> > &StartValueVec,
	  const vector<float> &StartDeltas, uint TriesPerIter,
	  uint MaxIters, uint MaxFailIters, float Shrink);
	void GetCloud(uint Index, uint MaxCloudSize,
	  float MaxCloudScoreDist, float MaxParamDist, 
	  vector<uint> &Indexes) const;
	bool SpatterIter();
	void SpatterUpdateDeltas_Shrink(float Shrink);
	void LogSpatterDeltas() const;
	void PrintfSpatterDeltas() const;
	float GetScore(uint Index) const;
	float GetQ(uint Index) const;
	float GetTC(uint Index) const;
	const vector<float> &GetParams(uint Index) const;
	float GetParamDist(uint Index1, uint Index2) const;
	float GetScoreDiff(uint Index1, uint Index2) const;
	float GetScoreDist(uint Index1, uint Index2) const;
	};
