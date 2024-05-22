#pragma once

#include <random>

class M3AlnParams
	{
public:
	float m_SubstMx_Letter[20][20];

public:
	float m_GapOpen = FLT_MAX;
	bool m_CenterAdded = false;

	uint m_PerturbSeed = 0;
	float m_PerturbSubstMxDelta = 0;
	float m_PerturbGapParamsDelta = 0;
	float m_PerturbDistMxDelta = 0;

	bool m_PerturbSubstMxDone = false;
	bool m_PerturbGapParamsDone = false;

	float m_NucMatchScore = FLT_MAX;
	float m_NucMismatchScore = FLT_MAX;
	float m_TermGapOpen = FLT_MAX;
	float m_TermGapExt = FLT_MAX;
	bool m_Ready = false;
	string m_Linkage = "min";
	uint m_TreeIters = 1;
	string m_KmerDist = "66";

private:
	float m_Center = FLT_MAX;
	mutable minstd_rand m_MinStdRand;

public:
	void LogMe() const { Print(g_fLog); }
	void Print(FILE *f) const;

	void SetFromCmdLine(bool IsNucleo, bool DoLog = true);

	void SetBlosum(uint PctId, uint n,
	  float GapOpen, float Center,
	  uint PerturbSeed, float PerturbSubstMxDelta,
	  float PerturbGapParamsDelta, float PerturbDistMxDelta,
	  bool DoLog = true);

	void UpdateMx(const Mx2020 &Mx, float GapOpen, float Center,
	  bool DoLog);

	void PerturbMyParams();
	void PerturbSubstMx();
	void PerturbGapParams();
	void PerturbDistMx(vector<vector<float> > &DistMx) const;

private:
	void Perturb1(float &Param, float MaxDelta);
	void Perturb1(float &Param, float MaxDelta) const;
	void AddCenter(float x);
	uint GetRand() const;
	void InitPerturb(uint Seed);
	};
