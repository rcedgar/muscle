#pragma once

#include "multisequence.h"
#include <unordered_map>

class Mega
	{
public:
	static string m_FileName;
	static vector<string> m_Lines;
	static vector<string> m_FeatureNames;
	static vector<float> m_Weights;
	static vector<uint> m_AlphaSizes;
	static unordered_map<string, uint> m_LabelToIdx;
	static unordered_map<string, uint> m_SeqToIdx;

// log(P_i) for each letter (for HMM Insert states)
	static vector<vector<float> > m_LogProbsVec;

// log(P_ij) for each letter pair (for HMM Match state)
	static vector<vector<vector<float> > > m_LogProbMxVec;

// log-odds score matrix for S-W / N-W
	static vector<vector<vector<float> > > m_LogOddsMxVec;

	static vector<string> m_Labels;
	static vector<vector<vector<byte> > > m_Profiles;
	static vector<string> m_Seqs;
	static uint m_NextLineNr;
	static uint m_FeatureCount;
	static bool m_Loaded;
	static float m_GapOpen;
	static float m_GapExt;

public:
	static void FromMSA_AAOnly(const MultiSequence &Aln,
	  float GapOpen, float GapExt);
	static void FromFile(const string &FileName);
	static uint GetProfileCount() { return SIZE(m_Profiles); }
	static const vector<vector<byte> > &GetProfile(uint ProfileIdx);
	static const string &GetLabel(uint ProfileIdx);
	static uint GetFeatureCount() { return m_FeatureCount; }
	static uint GetAlphaSize(uint FeatureIndex);
	static float GetWeight(uint FeatureIndex);
	static const string &GetFeatureName(uint FeatureIndex);
	static const string &GetNextLine();
	static void GetNextFields(vector<string> &Fields,
	  uint ExpectedNrFields = UINT_MAX);
	static float GetInsScore(const vector<vector<byte> > &Profile, uint Pos);
	static float GetMatchScore(
	  const vector<vector<byte> > &ProfileX, uint PosX,
	  const vector<vector<byte> > &ProfileY, uint PosY);
	static float GetMatchScore_LogOdds(
	  const vector<vector<byte> > &ProfileX, uint PosX,
	  const vector<vector<byte> > &ProfileY, uint PosY);
	static void CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
	  vector<float> &Freqs);
	static void LogFeatureParams(uint Idx);
	static void LogMx(const string &Name, const vector<vector<float> > &Mx);
	static void LogVec(const string &Name, const vector<float> &Vec);
	static void AssertSymmetrical(const vector<vector<float> > &Mx);
	static void CalcFwdFlat_mega(
	  const vector<vector<byte> > &ProfileX,
	  const vector<vector<byte> > &ProfileY, float *Flat);
	static void CalcBwdFlat_mega(
	  const vector<vector<byte> > &ProfileX,
	  const vector<vector<byte> > &ProfileY, float *Flat);
	static uint GetAAFeatureIdx();

public:
	static uint GetGSIByLabel(const string &Label);
	static const string &GetLabelByGSI(uint GSI);
	static const vector<vector<byte> > *GetProfileByGSI(uint GSI);
	static const vector<vector<byte> > *GetProfileByLabel(const string &Label);
	static const vector<vector<byte> > *GetProfileBySeq(const string &Seq,
	  bool FailOnError);
	};
