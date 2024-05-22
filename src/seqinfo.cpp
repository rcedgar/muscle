#include "muscle.h"
#include "seqinfo.h"
#include "alpha.h"
#include "objmgr.h"
//#include "fastq.h"

SeqInfo::SeqInfo()  : Obj(OT_SeqInfo)
	{
	m_Index = UINT_MAX;
	m_Label = 0;
	m_Seq = 0;
	m_LabelBuffer = 0;
	m_SeqBuffer = 0;
	m_L = 0;
	m_RevComp = false;
	m_LabelBytes = 0;
	m_MaxL = 0;
	m_MaxLabelBytes = 0;
	m_ORFNucSeq = 0;
	m_IsORF = false;
	m_ORFNucLo = UINT_MAX;
	m_ORFNucHi = UINT_MAX;
	m_ORFNucL = UINT_MAX;
	m_ORFFrame = 0;
	}

SeqInfo::~SeqInfo()
	{
	if (m_SeqBuffer != 0)
		myfree(m_SeqBuffer);
	if (m_LabelBuffer != 0)
		myfree(m_LabelBuffer);
	}

void SeqInfo::OnZeroRefCount()
	{
	m_Index = UINT_MAX;
	m_Seq = 0;
	m_Label = 0;
	m_L = 0;
	m_RevComp = false;
	m_IsORF = false;
	}

void SeqInfo::Copy(const SeqInfo &rhs)
	{
	AllocSeq(rhs.m_L);

	m_Index = rhs.m_Index;
	m_L = rhs.m_L;

	unsigned LabelBytes = (unsigned) strlen(rhs.m_Label) + 1;
	AllocLabel(LabelBytes);

	memcpy(m_LabelBuffer, rhs.m_Label, LabelBytes);
	memcpy(m_SeqBuffer, rhs.m_Seq, rhs.m_L);

	m_Seq = m_SeqBuffer;
	m_Label = m_LabelBuffer;

	m_IsORF = rhs.m_IsORF;
	m_ORFFrame = rhs.m_ORFFrame;
	m_ORFNucLo = rhs.m_ORFNucLo;
	m_ORFNucHi = rhs.m_ORFNucHi;
	m_ORFNucL = rhs.m_ORFNucL;
	m_ORFNucSeq = rhs.m_ORFNucSeq;
	if (rhs.m_ORFNucSeq != 0)
		ObjMgr::Up(m_ORFNucSeq);
	}

void SeqInfo::Init(unsigned Index)
	{
	m_Index = Index;
	m_L = 0;
	m_LabelBytes = 0;
	m_Seq = m_SeqBuffer;
	m_Label = m_LabelBuffer;
	m_IsORF = false;
	m_RevComp = false;
	}

void SeqInfo::AllocLabel(unsigned n)
	{
	if (n <= m_MaxLabelBytes)
		return;

	unsigned NewMaxLabelBytes = n + 128;
	char *NewLabelBuffer = myalloc(char, NewMaxLabelBytes);
	myfree(m_LabelBuffer);
	m_LabelBuffer = NewLabelBuffer;
	m_Label = NewLabelBuffer;
	m_MaxLabelBytes = NewMaxLabelBytes;
	}

void SeqInfo::AllocSeq(unsigned n)
	{
	if (n < m_MaxL)
		{
		m_Seq = m_SeqBuffer;
		return;
		}

//	StartTimer(SI_Realloc);
	unsigned NewMaxL = UINT_MAX;
	if (n < 10000)
		NewMaxL = RoundUp(2*n + 4096, 4096);
	else if (n < 1000000)
		NewMaxL = RoundUp(2*n + 65536, 65536);
	else
		NewMaxL = RoundUp((3*n)/2 + 1048576, 1048576);
	if (NewMaxL < m_L || NewMaxL < n)
		{
		Warning("SeqInfo::AllocSeq(n=%u), m_L=%u, NewMaxL=%u\n",
		  n, m_L, NewMaxL);
		NewMaxL = n + 4096;
		}
	byte *NewSeqBuffer = myalloc(byte, NewMaxL);
	if (m_L > 0)
		memcpy(NewSeqBuffer, m_Seq, m_L);
	myfree(m_SeqBuffer);
	m_Seq = NewSeqBuffer;
	m_SeqBuffer = NewSeqBuffer;
	m_MaxL = NewMaxL;
//	EndTimer(SI_Realloc);
	}

void SeqInfo::Pad(unsigned L, char c, char q)
	{
	if (L <= m_L)
		return;
	AllocSeq(L);

	for (unsigned i = m_L; i < L; ++i)
		{
		m_SeqBuffer[i] = c;
		}
	m_L = L;
	}

void SeqInfo::SetLabel(const char *Label)
	{
	unsigned n = (unsigned) strlen(Label) + 1;
	AllocLabel(n);
	m_LabelBytes = n;
	memcpy(m_LabelBuffer, Label, n);
	m_Label = m_LabelBuffer;
	}

void SeqInfo::SetCopy(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	SetLabel(Label);
	AllocSeq(L);
	memcpy(m_SeqBuffer, Seq, L);
	m_Seq = m_SeqBuffer;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::SetPtrs(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	m_Label = Label;
	m_Seq = Seq;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::GetReverse(SeqInfo *RevSI) const
	{
	RevSI->AllocSeq(m_L);
	RevSI->SetLabel(m_Label);
	RevSI->m_Index = UINT_MAX;

	byte *RevSeq = RevSI->m_SeqBuffer;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		RevSeq[m_L - i - 1] = c;
		}

	RevSI->m_L = m_L;
	RevSI->m_RevComp = !m_RevComp;
	RevSI->m_Index = m_Index;
	}


void SeqInfo::LogMe() const
	{
	Log("SeqInfo(%p) Seq %p, Buff %p, L %u, MaxL %u >%s\n",
	  this,
	  m_Seq,
	  m_SeqBuffer,
	  m_L,
	  m_MaxL,
	  m_Label);
	Log("%*.*s\n", m_L, m_L, m_Seq);
	}

unsigned SeqInfo::GetMemBytes() const
	{
	return m_MaxLabelBytes + m_MaxL;
	}

unsigned SeqInfo::GetIL() const
	{
	if (m_IsORF)
		return m_ORFNucL;
	return m_L;
	}

void SeqInfo::ToFasta(FILE *f, const char *Label) const
	{
	if (m_L == 0)
		return;

	SeqToFasta(f, m_Seq, m_L, Label);
	}

void SeqInfo::ToFastx(FILE *f, const char *Label) const
	{
	ToFasta(f, Label);
	}

void SeqInfo::TruncateLength(unsigned L)
	{
	if (m_L >= L)
		m_L = L;
	}

void SeqInfo::StripRight(unsigned n)
	{
	asserta(n < m_L);
	m_L -= n;
	}

void SeqInfo::StripLeft(unsigned n)
	{
	asserta(n < m_L);
	m_L -= n;
	asserta(m_Seq == m_SeqBuffer);
	for (unsigned i = 0; i < m_L; ++i)
		m_SeqBuffer[i] = m_SeqBuffer[i+n];
	}

void SeqInfo::StripGaps()
	{
	AllocSeq(m_L);
	unsigned NewL = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		if (c == '.' || c == '-')
			continue;
		m_SeqBuffer[NewL++] = c;
		}
	m_L = NewL;
	}

unsigned SeqInfo::GetWildcardCount(bool Nucleo) const
	{
	const byte *CharToLetter = (Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	const unsigned AlphaSize = (Nucleo ? 4 : 20);
	unsigned Count = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		unsigned Letter = CharToLetter[c];
		if (Letter >= AlphaSize)
			++Count;
		}
	return Count;
	}

unsigned SeqInfo::GetNCount() const
	{
	unsigned Count = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		if (c == 'N' || c == 'n')
			++Count;
		}
	return Count;
	}
