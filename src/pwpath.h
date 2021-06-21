#ifndef	PWPath_h
#define PWPath_h

/***
Each PWEdge in a PWPath specifies a column in a pair-wise (PW) alignment.
"Path" is by analogy with the path through an HMM.
Edge types are:

	'M'		LetterA + LetterB
	'D'		LetterA + GapB
	'I'		GapB + LetterA

The mnemomic is Match, Delete, Insert (with respect to A).
Here is a global alignment of sequences A and B.

	A:	AMQT-F
	B:	-M-TIF

The path for this example is:

	Edge	cType	uPrefixLengthA	uPrefixLengthB
	0		D		1				0
	1		M		2				1
	2		D		3				1
	3		M		4				2
	4		I		4				3
	5		M		5				4

Given the starting positions in each alignment (e.g., column zero for
a global alignment), the prefix length fields are redundant; they are
included only for convenience and as a sanity check, we are not trying
to optimize for speed or space here. We use prefix lengths rather than
column indexes because of the problem of representing the special case
of a gap in the first position.
***/

class Seq;
class MSA;
class SatchmoParams;
class PW;
class TextFile;
class PWScore;

class PWEdge
	{
public:
	char cType;
	unsigned uPrefixLengthA;
	unsigned uPrefixLengthB;

	bool Equal(const PWEdge &e) const
		{
		return uPrefixLengthA == e.uPrefixLengthA &&
		  uPrefixLengthB == e.uPrefixLengthB &&
		  cType == e.cType;
		}
	};

class PWPath
	{
// Disable compiler defaults
private:
	PWPath &operator=(const PWPath &rhs);
	PWPath(const PWPath &rhs);

public:
	PWPath();
	virtual ~PWPath();

public:
	void Clear();
	void FromStr(const char Str[]);
	void Copy(const PWPath &Path);
	void AppendEdge(const PWEdge &Edge);
	void AppendEdge(char cType, unsigned uPrefixLengthA, unsigned uPrefixLengthB);
	void PrependEdge(const PWEdge &Edge);
	unsigned GetEdgeCount() const { return m_uEdgeCount; }
	const PWEdge &GetEdge(unsigned uEdgeIndex) const;
	void Validate(const PWScore &PWS) const;
	void Validate() const;
	void LogMe() const;
	void FromFile(TextFile &File);
	void ToFile(TextFile &File) const;
	void FromMSAPair(const MSA &msaA, const MSA &msaB);
	void AssertEqual(const PWPath &Path) const;
	bool Equal(const PWPath &Path) const;
	unsigned GetMatchCount() const;
	unsigned GetDeleteCount() const;
	unsigned GetInsertCount() const;

private:
	void ExpandPath(unsigned uAdditionalEdgeCount);

private:
	unsigned m_uEdgeCount;
	unsigned m_uArraySize;
	PWEdge *m_Edges;
	};

#endif	// PWPath_h
