#ifndef SeqVect_h
#define SeqVect_h

#include <vector>
#include "seq.h"

typedef std::vector<Seq *> SeqVectBase;

class SeqVect : public SeqVectBase
	{
public:
	SeqVect() {}
	virtual ~SeqVect();

private:
// Not implemented; prevent use of copy c'tor and assignment.
	SeqVect(const SeqVect &);
	SeqVect &operator=(const SeqVect &);

public:
	void FromFile(TextFile &File)
		{
		FromFASTAFile(File);
		}

	void FromFASTAFile(TextFile &File);
	void FromFASTAFile(const string &FileName);

	void ToFASTAFile(TextFile &File) const;
	void ToFASTAFile(const string &FileName) const;

	void ToFile(TextFile &File) const
		{
		ToFASTAFile(File);
		}

	void PadToMSA(MSA &msa);
	void Copy(const SeqVect &rhs);
	void StripGaps();
	void StripGapsAndWhitespace();
	void ToUpper();
	void Clear();
	unsigned Length() const { return (unsigned) size(); }
	unsigned GetSeqCount() const { return (unsigned) size(); }
	void AppendSeq(const Seq &s);
	void AppendString(const string &Label, const string &s);
	bool FindName(const char *ptrName, unsigned *ptruIndex) const;
	void LogMe() const;
	const char *GetSeqName(unsigned uSeqIndex) const;
	void GetLabel(uint SeqIndex, string &Label) const;
	unsigned GetSeqId(unsigned uSeqIndex) const;
	unsigned GetSeqIdFromName(const char *Name) const;
	unsigned GetSeqLength(unsigned uSeqIndex) const;
	void SetSeqId(unsigned uSeqIndex, unsigned uId);
	Seq &GetSeq(unsigned uIndex);
	Seq &GetSeqById(unsigned uId);
	const Seq &GetSeq(unsigned uIndex) const;

	ALPHA GuessAlpha() const;
	void FixAlpha();

#ifndef	_WIN32
	reference at(size_type i) { return operator[](i); }
	const_reference at(size_type i) const { return operator[](i); }
#endif
	};

#endif	// SeqVect_h
