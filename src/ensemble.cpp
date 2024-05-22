#include "muscle.h"
#include "ensemble.h"
#include "qscorer.h"

#pragma warning(3: 4365)	// signed/unsigned conversion

static char ReadFirstChar(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	char c;
	ReadStdioFile(f, &c, 1);
	CloseStdioFile(f);
	return c;
	}

void Ensemble::SetDerived()
	{
	ToUpper();
	MapLabels();
	SortMSAs();
	SetUngappedSeqs();
	SetColToPosVec();
	SetColumns();
	}

void Ensemble::MapLabels()
	{
	asserta(!m_MSAs.empty());
	const MSA &M0 = *m_MSAs[0];
	const uint SeqCount = M0.GetSeqCount();
	M0.GetLabelToSeqIndex(m_Labels0, m_LabelToSeqIndex0);
	asserta(SIZE(m_Labels0) == SeqCount);
	}

void Ensemble::SortMSA(MSA &M)
	{
	const MSA &M0 = *m_MSAs[0];
	asserta(&M != &M0);

	const uint SeqCount = GetSeqCount();
	map<string, uint> LabelToSeqIndex2;
	vector<string> Labels2;
	M.GetLabelToSeqIndex(Labels2, LabelToSeqIndex2);

	char **szSeqsSorted = myalloc(char *, SeqCount);
	memset_zero(szSeqsSorted, SeqCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Labels2[SeqIndex];
		map<string, uint>::const_iterator p = m_LabelToSeqIndex0.find(Label);
		if (p == m_LabelToSeqIndex0.end())
			Die("SortMSA, different labels (%s)", Label.c_str());
		uint SeqIndex0 = p->second;
		asserta(szSeqsSorted[SeqIndex0] == 0);
		szSeqsSorted[SeqIndex0] = M.m_szSeqs[SeqIndex];
		}
	M.m_szNames = M0.m_szNames;
	M.m_szSeqs = szSeqsSorted;

	M.GetLabelToSeqIndex(Labels2, LabelToSeqIndex2);
	asserta(Labels2 == m_Labels0);
	asserta(LabelToSeqIndex2 == m_LabelToSeqIndex0);
	}

void Ensemble::SortMSAs()
	{
	const uint MSACount = GetMSACount();
	const uint SeqCount = GetSeqCount();

	const MSA &M0 = *m_MSAs[0];
	for (uint MSAIndex = 1; MSAIndex < MSACount; ++MSAIndex)
		{
		MSA &M = *m_MSAs[MSAIndex];
		const uint SeqCount2 = M.GetSeqCount();
		if (SeqCount2 != SeqCount)
			Die("Bad ensemble, different nr seqs");
		SortMSA(M);
		}
	}

void Ensemble::FromEFA(const string &FN)
	{
	Clear();

	vector<string> Strings;
	ReadStringsFromFile(FN, Strings);
	if (Strings.empty())
		Die("Empty EFA (%s)", FN.c_str());
	if (Strings[0].c_str()[0] != '<')
		Die("Invalid EFA, must start with '<' (%s)", FN.c_str());

	vector<string> MSAStrings;
	for (uint i = 0; i < SIZE(Strings); ++i)
		{
		const string &s = Strings[i];
		if (s.c_str()[0] == '<')
			{
			if (!MSAStrings.empty())
				{
				MSA &M = *new MSA;
				M.FromStrings(MSAStrings);
				m_MSAs.push_back(&M);
				MSAStrings.clear();
				}
			string MSAName = s.substr(1);
			m_MSANames.push_back(MSAName);
			}
		else
			MSAStrings.push_back(s);
		}

	MSA &M = *new MSA;
	M.FromStrings(MSAStrings);
	m_MSAs.push_back(&M);
	if (SIZE(m_MSAs) != SIZE(m_MSANames))
		Die("Invalid EFA, %u MSAs %u names (%s)",
		  SIZE(m_MSAs), SIZE(m_MSANames), FN.c_str());
	SetDerived();
	}

void Ensemble::ToEFA(const string &FN) const
	{
	if (FN.empty())
		return;
	FILE *f = CreateStdioFile(FN);
	const uint MSACount = GetMSACount();
	asserta(SIZE(m_MSANames) == MSACount);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &MSAName = m_MSANames[MSAIndex];
		fprintf(f, "<%s\n", MSAName.c_str());
		const MSA &M = *m_MSAs[MSAIndex];
		M.ToFASTAFile(f);
		}
	CloseStdioFile(f);
	}

void Ensemble::FromFile(const string &FN)
	{
	char c = ReadFirstChar(FN);
	if (c == '<')
		FromEFA(FN);
	else
		FromMSAPaths(FN);
	}

void Ensemble::FromMSAPaths(const string &FN)
	{
	Clear();

	m_MSANames.clear();
	ReadStringsFromFile(FN, m_MSANames);
	const uint MSACount = SIZE(m_MSANames);
	if (MSACount == 0)
		{
		Warning("Empty ensemble");
		return;
		}

	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		ProgressStep(MSAIndex, MSACount, "Reading m_MSAs");
		const string &MSAFileName = m_MSANames[MSAIndex];
		MSA *M = new MSA;
		M->FromFASTAFile(MSAFileName);
		m_MSAs.push_back(M);
		}

	if (opt(basename))
		{
		for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
			{
			const string &MSAFileName = m_MSANames[MSAIndex];
			m_MSANames[MSAIndex] = string(BaseName(MSAFileName.c_str()));
			}
		}

	if (opt(intsuffix))
		{
		for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
			{
			string Name = m_MSANames[MSAIndex];
			Psa(Name, ".%u", MSAIndex);
			m_MSANames[MSAIndex] = Name;
			}
		}

	SetDerived();
	}

uint Ensemble::GetSeqCount() const
	{
	if (m_MSAs.empty())
		return 0;
	uint SeqCount = m_MSAs[0]->GetSeqCount();
	return SeqCount;
	}

void Ensemble::ToUpper()
	{
	const uint MSACount = GetMSACount();
	const uint SeqCount = GetSeqCount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		MSA &M = *m_MSAs[MSAIndex];
		const uint ColCount = M.GetColCount();
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char *Seq = M.m_szSeqs[SeqIndex];
			for (uint i = 0; i < ColCount; ++i)
				Seq[i] = toupper(Seq[i]);
			}
		}
	}

void Ensemble::MakeResampledMSA(const vector<uint> &UniqueIxs, MSA &M) const
	{
	M.Clear();
	const uint ColCount = SIZE(UniqueIxs);
	const uint SeqCount = GetSeqCount();
	M.SetSize(SeqCount, ColCount);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = m_Labels0[SeqIndex];
		M.m_szNames[SeqIndex] = mystrsave(Label.c_str());
		}

	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint UniqueIx = UniqueIxs[ColIndex];
		asserta(UniqueIx < SIZE(m_UniqueIxs));
		uint Ix = m_UniqueIxs[UniqueIx];
		asserta(Ix < SIZE(m_ColumnStrings));
		const string &ColumnString = m_ColumnStrings[Ix];
		asserta(SIZE(ColumnString) == SeqCount);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char c = ColumnString[SeqIndex];
			M.m_szSeqs[SeqIndex][ColIndex] = c;
			}
		}
	}

void Ensemble::GetHiQualUniqueIxs(double MaxGapFract, double MinConf,
  vector<uint> &UniqueIxs) const
	{
	UniqueIxs.clear();
	const uint N = SIZE(m_UniqueIxs);
	for (uint UniqueIx = 0; UniqueIx < N; ++UniqueIx)
		{
		uint Ix = m_UniqueIxs[UniqueIx];
		double Conf = GetConf(UniqueIx);
		if (Conf < MinConf)
			continue;
		double GapFract = GetGapFract(Ix);
		if (GapFract <= MaxGapFract)
			UniqueIxs.push_back(UniqueIx);
		}
	}

uint Ensemble::GetMedianHiQualColCount(double MaxGapFract, double MinConf) const
	{
	vector<uint> ColCounts;
	const uint MSACount = SIZE(m_MSAs);
	if (MSACount == 0)
		return 0;

	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *m_MSAs[MSAIndex];
		const uint SeqCount = M.GetSeqCount();
		const uint ColCount = M.GetColCount();
		uint NonGappyColCount = 0;
		for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			double Conf = GetConf_MSACol(MSAIndex, ColIndex);
			if (Conf < MinConf)
				continue;
			uint GapCount = M.GetGapCount(ColIndex);
			double GapFract = double(GapCount)/double(SeqCount);
			if (GapFract <= MaxGapFract)
				++NonGappyColCount;
			}
		ColCounts.push_back(NonGappyColCount);
		}
	sort(ColCounts.begin(), ColCounts.end());
	uint MedianColCount = ColCounts[MSACount/2];
	return MedianColCount;
	}

void Ensemble::SetUngappedSeqs()
	{
	m_UngappedSeqs.clear();
	const MSA &M0 = *m_MSAs[0];
	const uint SeqCount = M0.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string UngappedSeq;
		M0.GetUngappedSeqStr(SeqIndex, UngappedSeq);
		m_UngappedSeqs.push_back(UngappedSeq);
		}

// Validate same seqs
	const uint MSACount = GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *m_MSAs[MSAIndex];
		asserta(M.GetSeqCount() == SeqCount);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			asserta(strcmp(M0.m_szNames[SeqIndex], M.m_szNames[SeqIndex]) == 0);
			string UngappedSeq;
			M.GetUngappedSeqStr(SeqIndex, UngappedSeq);
			if (UngappedSeq != m_UngappedSeqs[SeqIndex])
				{
				const uint L = SIZE(UngappedSeq);
				const uint L2 = SIZE(m_UngappedSeqs[SeqIndex]);
				Log(">%s\n", M0.m_szNames[SeqIndex]);
				Log("%s\n", UngappedSeq.c_str());
				Log("%s\n", m_UngappedSeqs[SeqIndex].c_str());
				for (uint i = 0; i < max(L, L2); ++i)
					{
					if (i >= min(L, L2))
						{
						Log("*");
						continue;
						}
					char c = UngappedSeq[i];
					char c2 = m_UngappedSeqs[SeqIndex][i];
					if (c == c2)
						Log(" ");
					else
						Log("d");
					}
				Log("\n");
				Die("MSA %u UngappedSeq != m_UngappedSeqs[%u]",
				  MSAIndex, SeqIndex);
				}
			}
		}
	}

void Ensemble::SetColToPosVec()
	{
	const uint MSACount = GetMSACount();
	const uint SeqCount = GetSeqCount();
	m_ColToPosVec.clear();
	m_ColToPosVec.resize(MSACount);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *m_MSAs[MSAIndex];
		m_ColToPosVec[MSAIndex].resize(SeqCount);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			M.GetColToPos1(SeqIndex, m_ColToPosVec[MSAIndex][SeqIndex]);
		}
	}

void Ensemble::GetColumn(uint MSAIndex, uint ColIndex,
  string &ColStr, vector<int> &PosVec) const
	{
	ColStr.clear();
	PosVec.clear();

	const MSA &M = *m_MSAs[MSAIndex];
	const uint SeqCount = GetSeqCount();
	ColStr.resize(SeqCount, '?');
	PosVec.resize(SeqCount, INT_MAX);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = M.GetChar(SeqIndex, ColIndex);
		ColStr[SeqIndex] = c;

		int Pos = m_ColToPosVec[MSAIndex][SeqIndex][ColIndex];
		PosVec[SeqIndex] = Pos;

		if (Pos > 0)
			{
			const string &UngappedSeq = m_UngappedSeqs[SeqIndex];
			asserta(uint(Pos) <= UngappedSeq.size());
			char c2 = UngappedSeq[uint(Pos)-1];
			asserta(c2 == c);
			}
		}
#if DEBUG
	{
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		asserta(ColStr[SeqIndex] != '?');
	}
#endif
	}

void Ensemble::SetColumns()
	{
	m_ColumnStrings.clear();
	m_ColumnPositions.clear();
	m_IxToMSAIndex.clear();
	m_IxToColIndex.clear();

	const uint MSACount = GetMSACount();
	const uint SeqCount = GetSeqCount();
	if (MSACount == 0)
		return;

	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		map<string, uint> LabelToSeqIndex2;
		const MSA &M = *m_MSAs[MSAIndex];
		uint SeqCount2 = M.GetSeqCount();
		asserta(SeqCount2 == SeqCount);

		const uint ColCount = M.GetColCount();
		for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			string ColStr;
			vector<int> PosVec;
			GetColumn(MSAIndex, ColIndex, ColStr, PosVec);

			m_ColumnStrings.push_back(ColStr);
			m_ColumnPositions.push_back(PosVec);
			m_IxToMSAIndex.push_back(MSAIndex);
			m_IxToColIndex.push_back(ColIndex);
			}
		}
	SetUniqueColMap();
	}

void Ensemble::SetUniqueColMap()
	{
	m_UniqueIxs.clear();
	m_UniqueIxToIxs.clear();
	m_IxToUniqueIx.clear();
	m_UniqueColMap.clear();
	m_MSAColToIx.clear();

	const uint MSACount = GetMSACount();
	m_MSAColToIx.resize(MSACount);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *m_MSAs[MSAIndex];
		uint ColCount = M.GetColCount();
		m_MSAColToIx[MSAIndex].resize(ColCount, UINT_MAX);
		}
	
	const vector<uint> Empty;
	const uint N = SIZE(m_ColumnPositions);
	for (uint Ix = 0; Ix < N; ++Ix)
		{
		uint MSAIndex = m_IxToMSAIndex[Ix];
		uint ColIndex = m_IxToColIndex[Ix];
		asserta(MSAIndex < MSACount);
		asserta(ColIndex < SIZE(m_MSAColToIx[MSAIndex]));
		m_MSAColToIx[MSAIndex][ColIndex] = Ix;
		const vector<int> &PosVec = m_ColumnPositions[Ix];
		map<vector<int>, uint>::const_iterator p =
		  m_UniqueColMap.find(PosVec);
		if (p == m_UniqueColMap.end())
			{
			uint UniqueIx = SIZE(m_UniqueIxs);
			m_UniqueColMap[PosVec] = UniqueIx;
			m_UniqueIxs.push_back(Ix);

			asserta(SIZE(m_UniqueIxToIxs) == UniqueIx);
			m_UniqueIxToIxs.push_back(Empty);
			m_UniqueIxToIxs[UniqueIx].push_back(Ix);
			m_IxToUniqueIx.push_back(UniqueIx);
			}
		else
			{
			uint UniqueIx = p->second;
			m_UniqueIxToIxs[UniqueIx].push_back(Ix);
			m_IxToUniqueIx.push_back(UniqueIx);
			}
		}
	ValidateUniqueColMap();
	}

void Ensemble::ValidateUniqueColMap1(uint MSAIndex, uint ColIndex) const
	{
	asserta(ColIndex < SIZE(m_MSAColToIx[MSAIndex]));
	uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
	asserta(Ix < SIZE(m_ColumnPositions));
	const vector<int> &PosVec = m_ColumnPositions[Ix];

	asserta(Ix < SIZE(m_IxToUniqueIx));
	uint UniqueIx = m_IxToUniqueIx[Ix];

	asserta(UniqueIx < SIZE(m_UniqueIxToIxs));
	const vector<uint> &Ixs = m_UniqueIxToIxs[UniqueIx];
	bool Found = false;
	for (uint i = 0; i < SIZE(Ixs); ++i)
		{
		if (Ixs[i] == Ix)
			{
			Found = true;
			break;
			}
		}
	asserta(Found);

	map<vector<int>, uint>::const_iterator p =
	  m_UniqueColMap.find(PosVec);
	asserta(p != m_UniqueColMap.end());
	uint UniqueIx2 = p->second;
	asserta(UniqueIx == UniqueIx2);
	}

void Ensemble::ValidateUniqueIx(uint UniqueIx) const
	{
	asserta(UniqueIx < SIZE(m_UniqueIxs));
	asserta(UniqueIx < SIZE(m_UniqueIxToIxs));

	uint Ix = m_UniqueIxs[UniqueIx];
	asserta(Ix < SIZE(m_ColumnPositions));

	const vector<int> &PosVec = m_ColumnPositions[Ix];
	map<vector<int>, uint>::const_iterator p =
	  m_UniqueColMap.find(PosVec);
	asserta(p != m_UniqueColMap.end());
	asserta(p->first == PosVec);

	const vector<uint> &Ixs = m_UniqueIxToIxs[UniqueIx];
	for (uint i = 0; i < SIZE(Ixs); ++i)
		{
		uint Ix2 = Ixs[i];
		asserta(Ix2 < SIZE(m_IxToUniqueIx));
		uint UniqueIx2 = m_IxToUniqueIx[Ix2];
		asserta(UniqueIx2 == UniqueIx);

		const vector<int> &PosVec2 = m_ColumnPositions[Ix2];
		asserta(PosVec2 == PosVec);
		}
	}

void Ensemble::ValidateUniqueColMap() const
	{
	const uint MSACount = GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = *m_MSAs[MSAIndex];
		uint ColCount = M.GetColCount();
		asserta(MSAIndex < SIZE(m_MSAColToIx));
		for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			ValidateUniqueColMap1(MSAIndex, ColIndex);
		}

	const uint UniqueIxCount = SIZE(m_UniqueIxs);
	for (uint UniqueIx = 0; UniqueIx < UniqueIxCount; ++UniqueIx)
		ValidateUniqueIx(UniqueIx);
	}

double Ensemble::GetGapFract(uint Ix) const
	{
	asserta(Ix < SIZE(m_ColumnStrings));
	const string &ColStr = m_ColumnStrings[Ix];
	const uint SeqCount = GetSeqCount();
	asserta(SIZE(ColStr) == SeqCount);
	uint GapCount = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		char c = ColStr[i];
		if (isgap(c))
			++GapCount;
		}
	double GapFract = double(GapCount)/SeqCount;
	return GapFract;
	}

void Ensemble::SubsampleWithReplacement(double MaxGapFract,
  uint ColCount, MSA &M) const
	{
	vector<uint> Ixs;
	GetIxSubset(MaxGapFract, Ixs);
	SubsampleWithReplacement(Ixs, ColCount, M);
	}

void Ensemble::SubsampleWithReplacement(const vector<uint> &Ixs,
  uint ColCount, MSA &M) const
	{
	asserta(ColCount > 0);
	const uint SeqCount = GetSeqCount();
	M.SetSize(SeqCount, ColCount);

	asserta(SIZE(m_Labels0) == SeqCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = m_Labels0[SeqIndex];
		M.m_szNames[SeqIndex] = mystrsave(Label.c_str());
		}

	const uint N = SIZE(Ixs);
	for (uint i = 0; i < N; ++i)
		{
		uint r = randu32()%N;
		uint Ix = Ixs[r];
		asserta(Ix < SIZE(m_ColumnStrings));
		const string &ColStr = m_ColumnStrings[Ix];
		asserta(SIZE(ColStr) == SeqCount);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			M.m_szSeqs[SeqIndex][i] = ColStr[SeqIndex];
		}
	}

void Ensemble::GetIxSubset(double MaxGapFract, vector<uint> &Ixs) const
	{
	Ixs.clear();
	const uint IxCount = SIZE(m_ColumnStrings);
	for (uint Ix = 0; Ix < IxCount; ++Ix)
		{
		double GapFract = GetGapFract(Ix);
		if (GapFract <= MaxGapFract)
			Ixs.push_back(Ix);
		}
	}


void Ensemble::GetAbToCountAll(vector<uint> &AbToCountAll)
	{
	const uint MSACount = GetMSACount();
	AbToCountAll.clear();
	AbToCountAll.resize(MSACount+1, 0);
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		vector<uint> AbToCount;
		GetAbToCount(MSAIndex, AbToCount);
		asserta(SIZE(AbToCount) == MSACount);
		for (uint i = 0; i < MSACount; ++i)
			AbToCountAll[i] += AbToCount[i];
		}
	}

void Ensemble::GetAbToCount(uint MSAIndex, vector<uint> &AbToCount)
	{
	const uint MSACount = GetMSACount();
	asserta(MSAIndex < MSACount);
	AbToCount.clear();
	AbToCount.resize(MSACount+1, 0);
	const MSA &M = *m_MSAs[MSAIndex];
	const uint ColCount = M.GetColCount();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint Ab = GetAb(MSAIndex, Col);
		asserta(Ab > 0);
		asserta(Ab <= MSACount);
		++AbToCount[Ab];
		}
	}

uint Ensemble::GetIx(uint MSAIndex, uint ColIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAColToIx));
	asserta(ColIndex < SIZE(m_MSAColToIx[MSAIndex]));
	uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
	return Ix;
	}

uint Ensemble::GetUniqueIx(uint MSAIndex, uint ColIndex) const
	{
	uint Ix = GetIx(MSAIndex, ColIndex);
	asserta(Ix < SIZE(m_IxToUniqueIx));
	uint UniqueIx = m_IxToUniqueIx[Ix];
	return UniqueIx;
	}

double Ensemble::GetMedianConf(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAs));
	const MSA &M = *m_MSAs[MSAIndex];
	const uint ColCount = M.GetColCount();
	vector<double> Confs;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
		asserta(Ix < SIZE(m_IxToUniqueIx));
		uint UniqueIx = m_IxToUniqueIx[Ix];
		double Conf = GetConf(UniqueIx);
		Confs.push_back(Conf);
		}
	sort(Confs.begin(), Confs.end());
	double MedianConf = Confs[ColCount/2];
	return MedianConf;
	}

double Ensemble::GetTotalConf(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAs));
	const MSA &M = *m_MSAs[MSAIndex];
	const uint ColCount = M.GetColCount();
	asserta(MSAIndex < SIZE(m_MSAColToIx));
	asserta(SIZE(m_MSAColToIx[MSAIndex]) == ColCount);
	double SumConf = 0;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
		asserta(Ix < SIZE(m_IxToUniqueIx));
		uint UniqueIx = m_IxToUniqueIx[Ix];
		double Conf = GetConf(UniqueIx);
		SumConf += Conf;
		}
	return SumConf;
	}

double Ensemble::GetConf_MSACol(uint MSAIndex, uint ColIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAColToIx));
	asserta(ColIndex < SIZE(m_MSAColToIx[MSAIndex]));
	uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
	uint UniqueIx = m_IxToUniqueIx[Ix];
	double Conf = GetConf(UniqueIx);
	return Conf;
	}

double Ensemble::GetConf(uint UniqueIx) const
	{
	const uint MSACount = GetMSACount();
	asserta(UniqueIx < SIZE(m_UniqueIxToIxs));
	const vector<uint> &Ixs = m_UniqueIxToIxs[UniqueIx];
	uint Ab = SIZE(Ixs);
	double Conf = double(Ab)/MSACount;
	return Conf;
	}

uint Ensemble::GetAb(uint MSAIndex, uint ColIndex) const
	{
	uint UniqueIx = GetUniqueIx(MSAIndex, ColIndex);
	asserta(UniqueIx < SIZE(m_UniqueIxToIxs));
	const vector<uint> &Ixs = m_UniqueIxToIxs[UniqueIx];
	uint Ab = SIZE(Ixs);
	asserta(Ab > 0);
	return Ab;
	}

uint Ensemble::GetN1(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAs));
	const MSA &M = *m_MSAs[MSAIndex];
	const uint ColCount = M.GetColCount();
	uint N1 = 0;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint Ab = GetAb(MSAIndex, ColIndex);
		if (Ab == 1)
			++N1;
		}
	return N1;
	}

void Ensemble::GetDispersion(double MaxGapFract,
  double &D_LetterPairs, double &D_Columns) const
	{
	QScorer QS;
	QS.m_MaxGapFract = MaxGapFract;

	vector<double> Qs;
	vector<double> TCs;
	const uint MSACount = GetMSACount();
	const uint PairCount = (MSACount*(MSACount - 1))/2;
	uint PairIndex = 0;
	for (uint i = 0; i < MSACount; ++i)
		{
		const MSA &MSAi = *m_MSAs[i];
		const string &Namei = m_MSANames[i];
		for (uint j = i + 1; j < MSACount; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Pairwise dists");

			const MSA &MSAj = *m_MSAs[j];
			const string &Namej = m_MSANames[j];
			QS.Run(Namei, MSAi, MSAj);
			double Qij = QS.m_Q;
			double TCij = QS.m_TC;

			QS.Run(Namej, MSAj, MSAi);
			double Qji = QS.m_Q;
			double TCji = QS.m_TC;

			double Q = (Qij + Qji)/2;
			double TC = (TCij + TCji)/2;
			asserta(Q >= 0 && Q <= 1);
			asserta(TC >= 0 && TC <= 1);
			Qs.push_back(Q);
			TCs.push_back(TC);
			}
		}
	sort(Qs.begin(), Qs.end());
	sort(TCs.begin(), TCs.end());
	const uint N = SIZE(Qs);
	asserta(SIZE(TCs) == N);

	double MedianQ = Qs[N/2];
	double MedianTC = TCs[N/2];

	D_LetterPairs = 1.0 - MedianQ;
	D_Columns = 1.0 - MedianTC;
	asserta(D_LetterPairs >= 0 && D_LetterPairs <= 1);
	asserta(D_Columns >= 0 && D_Columns <= 1);
	}

void Ensemble::CheckRefMSA(const MSA &Ref) const
	{
	const uint SeqCount = Ref.GetSeqCount();
	const uint RefSeqCount = Ref.GetSeqCount();
	const uint RefColCount = Ref.GetColCount();

	if (RefSeqCount != SeqCount)
		Die("Different nr seqs");

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Label = string(Ref.GetSeqName(SeqIndex));
		if (Label != m_Labels0[SeqIndex])
			Die("GetRefUniqueIxs, not sorted");
		}
	}

void Ensemble::GetRefUniqueIxs(const MSA &Ref,
  set<uint> &UniqueIxs, double MaxGapFract) const
	{
	UniqueIxs.clear();
	CheckRefMSA(Ref);

	const uint SeqCount = GetSeqCount();
	const uint RefColCount = Ref.GetColCount();

	vector<vector<int> > ColToPosVec(SeqCount);
	for (uint RefSeqIndex = 0; RefSeqIndex < SeqCount; ++RefSeqIndex)
		Ref.GetColToPos1(RefSeqIndex, ColToPosVec[RefSeqIndex]);

	for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
		{
		bool IsUpper = Ref.ColIsUpper(RefColIndex, MaxGapFract);
		if (!IsUpper)
			continue;

		vector<int> PosVec(SeqCount, UINT_MAX);
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			asserta(SeqIndex < SIZE(ColToPosVec));
			asserta(RefColIndex < SIZE(ColToPosVec[SeqIndex]));
			int Pos = ColToPosVec[SeqIndex][RefColIndex];
			PosVec[SeqIndex] = Pos;
			}
		map<vector<int>, uint >::const_iterator p =
		  m_UniqueColMap.find(PosVec);
		if (p != m_UniqueColMap.end())
			{
			uint UniqueIx = p->second;
			UniqueIxs.insert(UniqueIx);
			}
		}
	}

void Ensemble::GetRefPosSet(const MSA &Ref, double MaxGapFract,
  set<pair<uint, int> > &PosSet) const
	{
	PosSet.clear();

	CheckRefMSA(Ref);
	const uint SeqCount = GetSeqCount();
	const uint RefColCount = Ref.GetColCount();

	vector<vector<int> > ColToPosVec(SeqCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		Ref.GetColToPos1(SeqIndex, ColToPosVec[SeqIndex]);

	for (uint RefColIndex = 0; RefColIndex < RefColCount; ++RefColIndex)
		{
		bool IsUpper = Ref.ColIsUpper(RefColIndex, MaxGapFract);
		if (!IsUpper)
			continue;
		
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			int Pos = ColToPosVec[SeqIndex][RefColIndex];
			pair<uint, int> SeqPos(SeqIndex, Pos);
			PosSet.insert(SeqPos);
			}
		}
	}

void Ensemble::GetTestUniqueIxs(uint MSAIndex,
  const set<pair<uint, int> > &RefPosSet, vector<uint> &UniqueIxs,
  vector<double> &Confs) const
	{
	UniqueIxs.clear();
	Confs.clear();

	const uint MSACount = GetMSACount();
	const uint SeqCount = GetSeqCount();
	asserta(MSAIndex < MSACount);

	const MSA &M = *m_MSAs[MSAIndex];
	const uint ColCount = M.GetColCount();
	asserta(MSAIndex < SIZE(m_MSAColToIx));
	asserta(SIZE(m_MSAColToIx[MSAIndex]) == ColCount);
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint Ix = m_MSAColToIx[MSAIndex][ColIndex];
		asserta(Ix < SIZE(m_ColumnPositions));
		const vector<int> &PosVec = m_ColumnPositions[Ix];
		asserta(SIZE(PosVec) == SeqCount);
		uint FoundCount = 0;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			int Pos = PosVec[SeqIndex];
			if (Pos <= 0)
				continue;

			pair<uint, int> SeqPos(SeqIndex, Pos);
			if (RefPosSet.find(SeqPos) != RefPosSet.end())
				++FoundCount;
			}
		if (FoundCount >= SeqCount/2)
			{
			asserta(Ix < SIZE(m_IxToUniqueIx));
			uint UniqueIx = m_IxToUniqueIx[Ix];
			double Conf = GetConf(UniqueIx);

			UniqueIxs.push_back(UniqueIx);
			Confs.push_back(Conf);
			}
		}
	}

const MSA &Ensemble::GetMSA(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSAs));
	return *m_MSAs[MSAIndex];
	}

const string &Ensemble::GetMSAName(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_MSANames));
	return m_MSANames[MSAIndex];
	}
