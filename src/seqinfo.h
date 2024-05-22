#ifndef seqinfo_h
#define seqinfo_h

#include "obj.h"
#include "gobuff.h"

class SeqInfo : public Obj
	{
	friend class ObjMgr;

public:
	unsigned m_Index;
	const char *m_Label;
	const byte *m_Seq;

// Buffers are non-zero iff memory owned by this object.
	char *m_LabelBuffer;
	byte *m_SeqBuffer;

	unsigned m_L;
	unsigned m_LabelBytes;
	unsigned m_MaxL;
	unsigned m_MaxLabelBytes;
	bool m_RevComp;
	bool m_IsORF;
	SeqInfo *m_ORFNucSeq;
	int m_ORFFrame;
	unsigned m_ORFNucLo;
	unsigned m_ORFNucHi;
	unsigned m_ORFNucL;

protected:
	SeqInfo();
	~SeqInfo();

public:
	virtual unsigned GetMemBytes() const;
	virtual void OnZeroRefCount();

public:
	unsigned GetIL() const;
	void ToFasta(FILE *f, const char *Label = 0) const;
	void ToFastq(FILE *f, const char *Label = 0) const;
	void ToFastx(FILE *f, const char *Label = 0) const;
	void LogMe() const;
	void Init(unsigned Index);
	void AllocLabel(unsigned n);
	void AllocSeq(unsigned n);
	void AllocQual(unsigned n);
	void SetCopy(unsigned Index, const char *Label, const byte *Seq, unsigned L);
	void SetPtrs(unsigned Index, const char *Label, const byte *Seq, unsigned L);
	void Copy(const SeqInfo &rhs);
	void SetLabel(const char *Label);
	void GetRevComp(SeqInfo *RC) const;
	void GetReverse(SeqInfo *RC) const;
	void RevCompInPlace();
	void TruncateQual(unsigned IntQual);
	void TruncateTail(unsigned IntQual);
	void TruncateLength(unsigned L);
	void StripLeft(unsigned n);
	void StripRight(unsigned n);
	void Pad(unsigned n, char c, char q = 0);
	void StripGaps();
	void ReplaceSize(unsigned Size);
	unsigned GetNCount() const;
	unsigned GetWildcardCount(bool Nucleo) const;
	char GetMinQualChar() const;
	byte GetMinIntQual() const;
	void GetQualStr(string &Qual) const;
	unsigned GetSize() const;
	};

#endif // seqinfo_h
