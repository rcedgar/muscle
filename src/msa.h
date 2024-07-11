#ifndef	MSA_h
#define MSA_h

struct PathEdge;
class TextFile;
class Seq;
class ClusterNode;
class NodeCounts;
class DataBuffer;
class Sequence;
class MultiSequence;

class MSA
	{
public:
	unsigned m_uSeqCount;
	unsigned m_uColCount;
	unsigned m_uCacheSeqLength;
	unsigned m_uCacheSeqCount;
	char **m_szSeqs;
	char **m_szNames;

	static unsigned m_uIdCount;

	unsigned *m_IdToSeqIndex;
	unsigned *m_SeqIndexToId;

public:
	MSA();
	virtual ~MSA();

public:
// Ways to create an MSA
	void FromStrings(const vector<string> &Strings);
	void FromStrings2(const vector<string> &Labels, vector<string> &Seqs);
	void FromFile(TextFile &File);
	void FromFASTAFile(const string &FileName);
	void FromFASTAFile_PreserveCase(const string &FileName);
	void FromFASTAFile(TextFile &File);
	void FromSeq(const Seq &s);
	void FromSequence(const Sequence &s);
	void FromMultiSequence(const MultiSequence &MS);
	void GetLabelToSeqIndex(vector<string> &Labels,
	  map<string, uint> &LabelToIndex) const;

	void ToFile(TextFile &File) const;
	void ToFASTAFile(TextFile &File) const;
	void ToFASTAFile(FILE *f) const;
	void ToFASTAFile(const string &FileName) const;

	void SetSize(unsigned uSeqCount, unsigned uColCount);
	void SetSeqCount(unsigned uSeqCount);
	char GetChar(unsigned uSeqIndex, unsigned uIndex) const;
	unsigned GetLetter(unsigned uSeqIndex, unsigned uIndex) const;
	//unsigned GetLetterEx(unsigned uSeqIndex, unsigned uIndex) const;
	const char *GetLabel(unsigned uSeqIndex) const { return GetSeqName(uSeqIndex); }
	const char *GetSeqName(unsigned uSeqIndex) const;
	void GetSeqLabel(uint SeqIndex, string &Label) const;
	unsigned GetSeqId(unsigned uSeqIndex) const;
	unsigned GetSeqIndex(unsigned uId) const;
	bool GetSeqIndex(unsigned uId, unsigned *ptruIndex) const;
	double GetOcc(unsigned uColIndex) const;
	bool IsGap(unsigned uSeqIndex, unsigned uColIndex) const;
	//bool IsWildcard(unsigned uSeqIndex, unsigned uColIndex) const;
	bool IsGapColumn(unsigned uColIndex) const;
	uint GetGapCount(uint ColIndex) const;
	void GetUpperLowerGapCount(uint ColIndex,
	  uint &NU, uint &NL, uint &NG, uint &NDots, uint &NDashes) const;
	bool ColumnHasGap(unsigned uColIndex) const;
	bool IsGapSeq(unsigned uSeqIndex) const;
	void GetUngappedSeqStr(uint SeqIndex, string &SeqStr) const;

	void SetChar(unsigned uSeqIndex, unsigned uColIndex, char c);
	void SetSeqName(unsigned uSeqIndex, const char szName[]);
	void SetSeqId(unsigned uSeqIndex, unsigned uId);
	bool HasGap() const;
	bool IsLegalLetter(unsigned uLetter) const;
	void GetSeq(unsigned uSeqIndex, Seq &seq) const;
	void GetRowStr(unsigned uSeqIndex, string &SeqStr) const;
	const char *GetSeqCharPtr(uint SeqIndex) const;
	void Copy(const MSA &msa);
	double GetCons(unsigned uColIndex) const;
	double GetAvgCons() const;
	double GetPctIdentityPair(unsigned uSeqIndex1, unsigned uSeqIndex2) const;
	double GetPctIdentityPair2(unsigned uSeqIndex1, unsigned uSeqIndex2) const;
	bool GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const;
	uint GetSeqIndex(const string &Label, bool FailOnError = true) const;
	void DeleteCol(unsigned uColIndex);
	void DeleteColumns(unsigned uColIndex, unsigned uColCount);
	void CopySeq(unsigned uToSeqIndex, const MSA &msaFrom, unsigned uFromSeqIndex);
	void DeleteSeq(unsigned uSeqIndex);
	bool IsEmptyCol(unsigned uColIndex) const;
	void DeleteAllGapCols(MSA &msa) const;

	//ALPHA GuessAlpha() const;
	//void FixAlpha();

	unsigned GetCharCount(unsigned uSeqIndex, unsigned uColIndex) const;
	const char *GetSeqBuffer(unsigned uSeqIndex) const;
	unsigned AlignedColIndexToColIndex(unsigned uAlignedColIndex) const;
	unsigned GetUngappedSeqLength(unsigned uSeqIndex) const;
	void GetPWID(unsigned uSeqIndex1, unsigned uSeqIndex2, double *ptrdPWID,
	  unsigned *ptruPosCount) const;

	void LogMe() const;

	double GetPctGroupIdentityPair(unsigned uSeqIndex1, unsigned uSeqIndex2) const;

	void Clear()
		{
		Free();
		}
	unsigned GetSeqCount() const
		{
		return m_uSeqCount;
		}
	unsigned GetColCount() const
		{
		return m_uColCount;
		}

	static bool SeqsEq(const MSA &a1, unsigned uSeqIndex1, const MSA &a2,
	  unsigned uSeqIndex2);
	void GetPosToCol(uint SeqIndex, vector<uint> &PosToCol) const;
	void GetColToPos(uint SeqIndex, vector<uint> &ColToPos) const;

// 1-based positions, if <0 the column has a gap in this
// sequence which opens at 1-based position (-Pos).
	void GetColToPos1(uint SeqIndex, vector<int> &ColToPos) const;
	bool ColIsUpper(uint ColIndex, double MaxGapFract) const;
	bool ColIsAligned(uint ColIndex) const;

	static void SetIdCount(unsigned uIdCount);

private:
	void Free();
	void AppendSeq(char *ptrSeq, unsigned uSeqLength, char *ptrLabel);
	void ExpandCache(unsigned uSeqCount, unsigned uColCount);
	void GetNameFromFASTAAnnotationLine(const char szLine[],
	  char szName[], unsigned uBytes);
	void CopyCol(unsigned uFromCol, unsigned uToCol);
	};

void DeleteGappedCols(MSA &msa);
void MSAFromColRange(const MSA &msaIn, unsigned uFromColIndex, unsigned uColCount,
  MSA &msaOut);
void MSACat(const MSA &msa1, const MSA &msa2, MSA &msaCat);
void MSAAppend(MSA &msa1, const MSA &msa2);
void MSAFromSeqSubset(const MSA &msaIn, const unsigned uSeqIndexes[], unsigned uSeqCount,
  MSA &msaOut);
void AssertMSAEq(const MSA &msa1, const MSA &msa2);
void AssertMSAEqIgnoreCaseAndGaps(const MSA &msa1, const MSA &msa2);
void MSASubsetByIds(const MSA &msaIn, const unsigned Ids[], unsigned uIdCount,
  MSA &msaOut);
void SetMSAWeightsMuscle(MSA &msa);
void SetClustalWWeightsMuscle(MSA &msa);
void SetThreeWayWeightsMuscle(MSA &msa);

#endif	// MSA_h
