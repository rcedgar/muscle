#include "muscle.h"
#include "omplock.h"

static MultiSequence *g_GlobalMS;
static uint g_GlobalMSSeqCount = 0;
static double g_GlobalMSMeanSeqLength = 0;
static uint g_GlobalMSMaxSeqLength = 0;
static unordered_map<string, uint> m_LabelToIdx;
static unordered_map<string, const Sequence *> m_LabelToSeq;

uint GetGSIByLabel(const string &Label)
	{
	unordered_map<string, uint>::const_iterator iter =
	  m_LabelToIdx.find(Label);
	if (iter == m_LabelToIdx.end())
		Die("GetGSIByLabel(%s)", Label.c_str());
	uint GSI = iter->second;
	return GSI;
	}

void GetLabelByGSI(uint GSI, string &Label)
	{
	Label = string(g_GlobalMS->GetLabel(GSI));
	}

uint GetSeqLengthByGSI(uint GSI)
	{
	uint L = g_GlobalMS->GetSeqLength(GSI);
	return L;
	}

uint GetSeqLengthByGlobalLabel(const string &Label)
	{
	//uint GSI = GetGSIByLabel(Label);
	//uint L = GetSeqLengthByGSI(GSI);
	const Sequence &seq = GetGlobalInputSeqByLabel(Label);
	uint L = seq.GetLength();
	return L;
	}

const Sequence *GetSequenceByGlobalLabel(const string &Label)
	{
	uint GSI = GetGSIByLabel(Label);
	const Sequence *seq = GetSequenceByGSI(GSI);
	return seq;
	}

const Sequence *GetSequenceByGSI(uint GSI)
	{
	return g_GlobalMS->GetSequence(GSI);
	}

const byte *GetByteSeqByGSI(uint GSI)
	{
	const Sequence *seq = GetSequenceByGSI(GSI);
	const byte *ByteSeq = seq->GetBytePtr();
	return ByteSeq;
	}

void AddGlobalTmpSeq(Sequence *seq)
	{
	Lock();
	string Label = string(seq->m_Label);
	m_LabelToSeq[Label] = seq;
	Unlock();
	}

void SetGlobalInputMS(MultiSequence &MS)
	{
	asserta(g_GlobalMS == 0);
	g_GlobalMS = &MS;
	g_GlobalMSSeqCount = g_GlobalMS->GetSeqCount();
	g_GlobalMSMeanSeqLength = 0;
	g_GlobalMSMaxSeqLength = 0;
	double SumSeqLength = 0;
	for (uint GSI = 0; GSI < g_GlobalMSSeqCount; ++GSI)
		{
		const Sequence *Seq = g_GlobalMS->GetSequence(GSI);
		string Label = string(g_GlobalMS->GetLabel(GSI));
		unordered_map<string, uint>::const_iterator iter =
		  m_LabelToIdx.find(Label);
		if (iter != m_LabelToIdx.end())
			Die("Error duplicate label in input >%s", Label.c_str());
		m_LabelToIdx[Label] = GSI;
		m_LabelToSeq[Label] = Seq;
		uint L = Seq->GetLength();
		g_GlobalMSMaxSeqLength = max(L, g_GlobalMSMaxSeqLength);
		SumSeqLength += L;
		}
	if (g_GlobalMSSeqCount > 0)
		g_GlobalMSMeanSeqLength = SumSeqLength/g_GlobalMSSeqCount;
	ShowSeqStats(*g_GlobalMS);
	}

MultiSequence &GetGlobalInputMS()
	{
	asserta(g_GlobalMS != 0);
	return *g_GlobalMS;
	}

uint GetGlobalMSSeqCount()
	{
	return g_GlobalMSSeqCount;
	}

double GetGlobalMSMeanSeqLength()
	{
	return g_GlobalMSMeanSeqLength;
	}

uint GetGSICount()
	{
	return GetGlobalMSSeqCount();
	}

const Sequence &GetGlobalInputSeqByIndex(uint GSI)
	{
	asserta(GSI < g_GlobalMSSeqCount);
	asserta(g_GlobalMS != 0);
	const Sequence *Seq = g_GlobalMS->GetSequence(GSI);
	asserta(Seq != 0);
	return *Seq;
	}

const Sequence &GetGlobalInputSeqByLabel(const string &Label)
	{
	//uint GSI = GetGSIByLabel(Label);
	//const Sequence &Seq = GetGlobalInputSeqByIndex(GSI);
	unordered_map<string, const Sequence *>::const_iterator iter = m_LabelToSeq.find(Label);
	if (iter == m_LabelToSeq.end())
		Die("GetGlobalInputSeqByLabel()", Label.c_str());
	const Sequence *seq = iter->second;
	const string &Label2 = seq->GetLabel();
	asserta(Label2 == Label);
	return *seq;
	}

const byte *GetGlobalByteSeqByLabel(const string &Label)
	{
	const Sequence &seq = GetGlobalInputSeqByLabel(Label);
	const byte *ByteSeq = seq.GetBytePtr();
	return ByteSeq;
	}
