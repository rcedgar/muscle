#pragma once

#include <set>

class Ensemble
	{
public:
	vector<MSA *> m_MSAs;
	vector<string> m_MSANames;
	vector<string> m_Labels0;
	map<string, uint> m_LabelToSeqIndex0;
	vector<string> m_UngappedSeqs;
	vector<string> m_ColumnStrings;
	vector<vector<int> > m_ColumnPositions;

// 1-based positions, if <0 the column has a gap in this
// sequence which opens at 1-based position (-Pos).
// If ColToPos[Col] is 0, this is left terminal gap.
	vector<vector<vector<int> > > m_ColToPosVec;

	vector<uint> m_IxToMSAIndex;
	vector<uint> m_IxToColIndex;
	vector<vector<uint> > m_MSAColToIx;

	vector<uint> m_UniqueIxs;
	vector<vector<uint> > m_UniqueIxToIxs;
	vector<uint> m_IxToUniqueIx;
	map<vector<int>, uint> m_UniqueColMap;

public:
	void Clear()
		{
		m_MSAs.clear();
		m_MSANames.clear();
		m_ColumnStrings.clear();
		m_ColumnPositions.clear();
		m_Labels0.clear();
		m_LabelToSeqIndex0.clear();
		m_ColToPosVec.clear();
		m_IxToMSAIndex.clear();
		m_IxToColIndex.clear();
		m_UniqueIxToIxs.clear();
		m_UniqueIxs.clear();
		m_IxToUniqueIx.clear();
		m_UniqueColMap.clear();
		m_MSAColToIx.clear();
		}

	void FromFile(const string &FileName);
	void FromMSAPaths(const string &FileName);
	void FromEFA(const string &FileName);

	void ToEFA(const string &FileName) const;

	void SetDerived();
	uint GetMSACount() const { return SIZE(m_MSAs); }
	uint GetIxCount() const { return SIZE(m_IxToMSAIndex); }
	uint GetSeqCount() const;
	void SetColumns();
	void GetColumn(uint MSAIndex, uint ColIndex,
	  string &ColStr, vector<int> &ColPos) const;
	void GetIxSubset(double MaxGapFract, vector<uint> &Ixs) const;
	double GetGapFract(uint Ix) const;
	void SubsampleWithReplacement(double MaxGapFract,
	  uint ColCount, MSA &M) const;
	void SubsampleWithReplacement(const vector<uint> &Ixs,
	  uint ColCount, MSA &M) const;
	void GetAbToCountAll(vector<uint> &AbToCount);
	void GetAbToCount(uint MSAIndex, vector<uint> &AbToCount);
	uint GetUniqueIx(uint MSAIndex, uint ColIndex) const;
	uint GetIx(uint MSAIndex, uint ColIndex) const;
	uint GetAb(uint MSAIndex, uint ColIndex) const;
	double GetConf(uint UniqueIx) const;
	double GetConf_MSACol(uint MSAIndex, uint ColIndex) const;
	uint GetN1(uint MSAIndex) const;
	void ValidateUniqueColMap() const;
	void ValidateUniqueColMap1(uint MSAIndex, uint ColIndex) const;
	void ValidateUniqueIx(uint UniqueIx) const;
	void GetDispersion(double MaxGapFract,
	  double &D_LetterPairs, double &D_Columns) const;
	double GetTotalConf(uint MSAIndex) const;
	double GetMedianConf(uint MSAIndex) const;

	void SortMSA(MSA &M);
	void CheckRefMSA(const MSA &Ref) const;
	void GetRefPosSet(const MSA &Ref, double MaxGapFract,
	  set<pair<uint, int> > &PosSet) const;
	void GetTestUniqueIxs(uint MSAIndex,
	  const set<pair<uint, int> > &RefPosSet, vector<uint> &UniqueIxs,
	  vector<double> &Confs) const;
	void GetRefUniqueIxs(const MSA &Ref, set<uint> &UniqueIxs,
	  double MaxGapFract) const;
	void MakeResampledMSA(const vector<uint> &UniqueIxs, MSA &M) const;

	uint GetMedianHiQualColCount(double MaxGapFract, double MinConf) const;
	void GetHiQualUniqueIxs(double MaxGapFract, double MinConf,
	  vector<uint> &UniqueIxs) const;

	const MSA &GetMSA(uint MSAIndex) const;
	const string &GetMSAName(uint MSAIndex) const;

	void GetLetterConfsVec(const MSA &Ref, double MaxGapFract,
	  vector<vector<double> > &LetterConfsVec) const;

private:
	void MapLabels();
	void SortMSAs();
	void ToUpper();
	void SetColToPosVec();
	void SetUngappedSeqs();
	void SetUniqueColMap();
	};
