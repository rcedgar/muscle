#include "myutils.h"
#include "mega.h"
#include "alpha.h"
#include "pairhmm.h"

void GetBlosum62LogOddsLetterMx(vector<vector<float> > &LogOddsMx);

string Mega::m_FileName;
vector<string> Mega::m_FeatureNames;
vector<float> Mega::m_Weights;
vector<uint> Mega::m_AlphaSizes;
vector<string> Mega::m_Labels;
vector<vector<vector<byte> > > Mega::m_Profiles;
vector<string> Mega::m_Seqs;
vector<vector<float> > Mega::m_LogProbsVec;
vector<vector<vector<float> > > Mega::m_LogOddsMxVec;
uint Mega::m_FeatureCount;
bool Mega::m_Loaded = false;
float Mega::m_GapOpen = FLT_MAX;
float Mega::m_GapExt = FLT_MAX;
unordered_map<string, uint> Mega::m_LabelToIdx;
unordered_map<string, uint> Mega::m_SeqToIdx;

void Mega::RejectLegacyMega(const string &FileName)
	{
	if (EndsWith(FileName, ".mega"))
		Die("Text .mega files are no longer supported; use STRUCTS "
		  "(.pdb/.cif/.cal/.bca/.bcb/.files or a directory)");
	}

uint Mega::GetGSIByLabel(const string &Label)
	{
	unordered_map<string, uint>::const_iterator iter = m_LabelToIdx.find(Label);
	if (iter == m_LabelToIdx.end())
		Die("Mega::GetGSIByLabel(%s)", Label.c_str());
	uint GSI = iter->second;
	return GSI;
	}

const string &Mega::GetLabelByGSI(uint GSI)
	{
	asserta(GSI < SIZE(m_Labels));
	return m_Labels[GSI];
	}

const vector<vector<byte> > *Mega::GetProfileByGSI(uint GSI)
	{
	asserta(GSI < SIZE(m_Profiles));
	return &m_Profiles[GSI];
	}

const vector<vector<byte> > *Mega::GetProfileByLabel(const string &Label)
	{
	unordered_map<string, uint>::const_iterator iter = m_LabelToIdx.find(Label);
	if (iter == m_LabelToIdx.end())
		Die("Mega::GetProfileByLabel(%s)", Label.c_str());
	uint Idx = iter->second;
	asserta(Idx < SIZE(m_Profiles));
	return &m_Profiles[Idx];
	}

const vector<vector<byte> > *Mega::GetProfileBySeq(const string &Seq,
  bool FailOnError)
	{
	unordered_map<string, uint>::const_iterator iter = m_SeqToIdx.find(Seq);
	if (iter == m_SeqToIdx.end())
		{
		if (FailOnError)
			Die("Mega::GetProfileBySeq(%16.16s...)", Seq.c_str());
		return 0;
		}
	uint Idx = iter->second;
	asserta(Idx < SIZE(m_Profiles));
	return &m_Profiles[Idx];
	}

void Mega::AssertSymmetrical(const vector<vector<float> > &Mx)
	{
	const uint N = SIZE(Mx);
	for (uint i = 0; i < N; ++i)
		{
		const vector<float> &Row = Mx[i];
		asserta(SIZE(Row) == N);
		for (uint j = 0; j < i; ++j)
			asserta(feq(Mx[i][j], Mx[j][i]));
		}
	}

void Mega::CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
  vector<float> &MarginalFreqs)
	{
	MarginalFreqs.clear();
	AssertSymmetrical(FreqsMx);
	const uint N = SIZE(FreqsMx);
	float SumMarginalFreqs = 0;
	for (uint i = 0; i < N; ++i)
		{
		const vector<float> &Row = FreqsMx[i];
		float MarginalFreq = 0;
		for (uint j = 0; j < N; ++j)
			MarginalFreq += Row[j];
		MarginalFreqs.push_back(MarginalFreq);
		SumMarginalFreqs += MarginalFreq;
		}
	asserta(feq(SumMarginalFreqs, 1));
	}

float Mega::GetInsScore(const vector<vector<byte> > &Profile, uint Pos)
	{
	asserta(Pos < SIZE(Profile));
	// STRUCTS / AA-only paths do not use insert emission scores.
	return 0.0f;
	}

const string &Mega::GetFeatureName(uint FeatureIndex)
	{
	asserta(FeatureIndex < SIZE(m_FeatureNames));
	return m_FeatureNames[FeatureIndex];
	}

uint Mega::GetAlphaSize(uint FeatureIndex)
	{
	asserta(FeatureIndex < SIZE(m_AlphaSizes));
	return m_AlphaSizes[FeatureIndex];
	}

float Mega::GetWeight(uint FeatureIndex)
	{
	asserta(FeatureIndex < SIZE(m_Weights));
	return m_Weights[FeatureIndex];
	}

const string &Mega::GetLabel(uint ProfileIdx)
	{
	asserta(ProfileIdx < SIZE(m_Profiles));
	return m_Labels[ProfileIdx];
	}

const vector<vector<byte> > &Mega::GetProfile(uint ProfileIdx)
	{
	asserta(ProfileIdx < SIZE(m_Profiles));
	return m_Profiles[ProfileIdx];
	}

float Mega::GetMatchScore_LogOdds(
  const vector<vector<byte> > &ProfileX, uint PosX,
  const vector<vector<byte> > &ProfileY, uint PosY)
	{
	const uint LX = SIZE(ProfileX);
	const uint LY = SIZE(ProfileY);
	asserta(PosX < LX);
	asserta(PosY < LY);
	const vector<byte> &ProfColX = ProfileX[PosX];
	const vector<byte> &ProfColY = ProfileY[PosY];
	float Score = 0;
	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		const vector<vector<float> > &SubstMx = m_LogOddsMxVec[i];
		byte LetterX = ProfColX[i];
		byte LetterY = ProfColY[i];
		float LetterPairScore = SubstMx[LetterX][LetterY];
		Score += LetterPairScore*m_Weights[i];
		}
	return Score;
	}

float Mega::GetMatchScore(
  const vector<vector<byte> > &ProfileX, uint PosX,
  const vector<vector<byte> > &ProfileY, uint PosY)
	{
	return GetMatchScore_LogOdds(ProfileX, PosX, ProfileY, PosY);
	}

void Mega::LogVec(const string &Name, const vector<float> &Vec)
	{
	const uint N = SIZE(Vec);
	Log("\n%s/%u", Name.c_str(), N);
	for (uint i = 0; i < N; ++i)
		{
		if (i%10 == 0)
			Log("\n  ");
		else
			Log(" ");
		Log("[%2u]=%.2f", i, Vec[i]);
		}
	Log("\n");
	}

void Mega::LogMx(const string &Name,
  const vector<vector<float> > &Mx)
	{
	const uint N = SIZE(Mx);
	Log("\n%s/%u\n", Name.c_str(), N);

	Log("     ");
	for (uint j = 0; j < N; ++j)
		Log(" %7u", j);
	Log("\n");
	for (uint i = 0; i < N; ++i)
		{
		Log("[%2u] ", i);
		const vector<float> &Row = Mx[i];
		asserta(SIZE(Row) == N);
		for (uint j = 0; j < N; ++j)
			Log(" %7.2f", Row[j]);
		Log("\n");
		}
	}

void Mega::LogFeatureParams(uint Idx)
	{
	asserta(Idx < SIZE(m_FeatureNames));
	asserta(Idx < SIZE(m_LogOddsMxVec));
	const string &Name = m_FeatureNames[Idx];
	Log("\n");
	Log("Feature %s, weight %.3g\n",
	  Name.c_str(), m_Weights[Idx]);
	if (Idx < SIZE(m_LogProbsVec) && !m_LogProbsVec[Idx].empty())
		LogVec(Name, m_LogProbsVec[Idx]);
	LogMx(Name, m_LogOddsMxVec[Idx]);
	}

bool Mega::IsAAFeatureName(const string &Name)
	{
	if (Name == "AA" || Name == "aa")
		return true;
	if (SIZE(Name) >= 3 &&
	  (Name[0] == 'A' || Name[0] == 'a') &&
	  (Name[1] == 'A' || Name[1] == 'a'))
		{
		for (uint i = 2; i < SIZE(Name); ++i)
			if (!isdigit(Name[i]))
				return false;
		return true;
		}
	return false;
	}

uint Mega::GetAAFeatureIdx()
	{
	for (uint FeatureIdx = 0; FeatureIdx < SIZE(m_FeatureNames); ++FeatureIdx)
		if (IsAAFeatureName(m_FeatureNames[FeatureIdx]))
			return FeatureIdx;
	Die("Mega::GetAAFeatureIdx(), not found");
	return UINT_MAX;
	}

void Mega::FromMSA_AAOnly(const MultiSequence &Aln,
  float GapOpen, float GapExt)
	{
	m_FileName = "FromMSA_AAOnly()";

	m_FeatureNames.clear();
	m_FeatureNames.push_back("AA");

	m_Weights.clear();
	m_Weights.push_back(1.0f);

	m_AlphaSizes.clear();
	m_AlphaSizes.push_back(20);
	m_FeatureCount = 1;

	m_LabelToIdx.clear();
	m_SeqToIdx.clear();
	m_Labels.clear();
	m_Seqs.clear();
	m_Profiles.clear();
	const uint SeqCount = Aln.GetSeqCount();
	m_Profiles.resize(SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		const string &Label = Aln.GetLabelStr(SeqIdx);
		m_Labels.push_back(Label);
		string Seq;
		Aln.GetSeqStr(SeqIdx, Seq);
		string UngappedSeq;
		for (uint i = 0; i < SIZE(Seq); ++i)
			{
			char c = Seq[i];
			if (!isgap(c))
				UngappedSeq += c;
			}
		m_Seqs.push_back(UngappedSeq);
		m_LabelToIdx[Label] = SeqIdx;
		m_SeqToIdx[UngappedSeq] = SeqIdx;

		vector<vector<byte> > &Profile = m_Profiles[SeqIdx];
		for (uint i = 0; i < SIZE(UngappedSeq); ++i)
			{
			char c = UngappedSeq[i];
			byte Letter = g_CharToLetterAmino[c];
			if (Letter >= 20)
				Letter = 0;
			vector<byte> Col;
			Col.push_back(Letter);
			Profile.push_back(Col);
			}
		}

	m_LogProbsVec.clear();
	m_LogOddsMxVec.clear();
	m_LogOddsMxVec.resize(1);
	GetBlosum62LogOddsLetterMx(m_LogOddsMxVec[0]);
	m_GapOpen = GapOpen;
	m_GapExt = GapExt;
	m_Loaded = true;
	}
