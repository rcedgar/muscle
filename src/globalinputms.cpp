#include "muscle.h"

static MultiSequence *g_GlobalMS;
static uint g_GlobalMSSeqCount = 0;
static double g_GlobalMSMeanSeqLength = 0;
static uint g_GlobalMSMaxSeqLength = 0;

void ClearGlobalInputMS()
	{
	if (g_GlobalMS == 0)
		return;
	delete g_GlobalMS;
	g_GlobalMS = 0;
	}

MultiSequence &LoadGlobalInputMS(const string &FileName)
	{
	asserta(g_GlobalMS == 0);
	g_GlobalMS = new MultiSequence;
	asserta(g_GlobalMS != 0);
	g_GlobalMS->FromFASTA(FileName, true);
	g_GlobalMSSeqCount = g_GlobalMS->GetSeqCount();
	g_GlobalMSMeanSeqLength = 0;
	g_GlobalMSMaxSeqLength = 0;
	double SumSeqLength = 0;
	for (uint GSI = 0; GSI < g_GlobalMSSeqCount; ++GSI)
		{
		const Sequence *Seq = g_GlobalMS->GetSequence(GSI);
		uint L = Seq->GetLength();
		g_GlobalMSMaxSeqLength = max(L, g_GlobalMSMaxSeqLength);
		SumSeqLength += L;
		Sequence *HackSeq = (Sequence *) Seq;
		HackSeq->m_GSI = GSI;
		}
	if (g_GlobalMSSeqCount > 0)
		g_GlobalMSMeanSeqLength = SumSeqLength/g_GlobalMSSeqCount;
	return *g_GlobalMS;
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

void ShowGlobalInputSeqStats()
	{
	ProgressLog("Input: %u seqs, length avg %.0f max %u\n\n",
	  g_GlobalMSSeqCount, g_GlobalMSMeanSeqLength, g_GlobalMSMaxSeqLength);

	if (g_GlobalMSMaxSeqLength > 5000)
		Warning("Sequence length >5k may require excessive memory");
	}
