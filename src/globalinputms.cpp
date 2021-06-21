#include "muscle.h"

static MultiSequence *g_GlobalMS;
static uint g_GlobalMSSeqCount = 0;

MultiSequence &LoadGlobalInputMS(const string &FileName)
	{
	asserta(g_GlobalMS == 0);
	g_GlobalMS = new MultiSequence;
	asserta(g_GlobalMS != 0);
	g_GlobalMS->FromFASTA(FileName);
	g_GlobalMSSeqCount = g_GlobalMS->GetSeqCount();
	for (uint GSI = 0; GSI < g_GlobalMSSeqCount; ++GSI)
		{
		Sequence *Seq = g_GlobalMS->GetSequence(GSI);
		Seq->SetGSI(GSI);
		}
	return *g_GlobalMS;
	}

const MultiSequence &GetGlobalInputMS()
	{
	asserta(g_GlobalMS != 0);
	return *g_GlobalMS;
	}

uint GetGlobalMSSeqCount()
	{
	return g_GlobalMSSeqCount;
	}

uint GetGSICount()
	{
	return GetGlobalMSSeqCount();
	}

const Sequence &GetGlobalInputSeq(uint GSI)
	{
	asserta(GSI < g_GlobalMSSeqCount);
	asserta(g_GlobalMS != 0);
	const Sequence *Seq = g_GlobalMS->GetSequence(GSI);
	asserta(Seq != 0);
	return *Seq;
	}

const string &GetGlobalInputSeqLabel(uint GSI)
	{
	const Sequence &Seq = GetGlobalInputSeq(GSI);
	const string &Label = Seq.GetLabel();
	return Label;
	}
