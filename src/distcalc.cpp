#include "muscle.h"
#include "distfunc.h"
#include "distcalc.h"
#include "msa.h"

double GetScoreDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2);
double GetMyDist(const MSA &msa, unsigned SeqIndex1, unsigned SeqIndex2);

void DistCalcDF::Init(const DistFunc &DF)
	{
	m_ptrDF = &DF;
	m_Type = DF.m_DistType;
	}

void DistCalcDF::CalcDistRange(unsigned i, dist_t Dist[]) const
	{
	for (unsigned j = 0; j < i; ++j)
		Dist[j] = m_ptrDF->GetDist(i, j);
	}

unsigned DistCalcDF::GetCount() const
	{
	return m_ptrDF->GetCount();
	}

unsigned DistCalcDF::GetId(unsigned i) const
	{
	return m_ptrDF->GetId(i);
	}

const char *DistCalcDF::GetName(unsigned i) const
	{
	return m_ptrDF->GetName(i);
	}

void DistCalcMSA::Init(const MSA &msa, DISTANCE Distance)
	{
	m_ptrMSA = &msa;
	m_Distance = Distance;
	m_Type = DISTANCEToStr(m_Distance);
	}

static float GetPctId(const char *Seq1, const char *Seq2, uint ColCount)
	{
	uint IdCount = 0;
	uint PosCount = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		const char c1 = Seq1[i];
		const char c2 = Seq2[i];
		if (IsGapChar(c1) || IsGapChar(c2))
			continue;
		if (c1 == c2)
			++IdCount;
		++PosCount;
		}
	if (0 == PosCount)
		return 0;
	return (float) IdCount / (float) PosCount;
	}

dist_t DistCalcMSA::CalcDist(uint i, uint j) const
	{
	dist_t d = 0;
	const uint ColCount = m_ptrMSA->GetColCount();
	const char *Seqi = m_ptrMSA->m_szSeqs[i];
	const char *Seqj = m_ptrMSA->m_szSeqs[j];
	switch (m_Distance)
		{
	case DISTANCE_PctIdKimura:
		{
		const float PctId = GetPctId(Seqi, Seqj, ColCount);
		d = (float) KimuraDist(PctId);
		break;
		}
	case DISTANCE_PctIdLog:
		{
		const float PctId = GetPctId(Seqi, Seqj, ColCount);
		d = (float) PctIdToMAFFTDist(PctId);
		break;
		}
	case DISTANCE_ScoreDist:
		{
		d = (float) GetScoreDist(*m_ptrMSA, i, j);
		break;
		}
	case DISTANCE_MyDist:
		{
		d = (float) GetMyDist(*m_ptrMSA, i, j);
		break;
		}
	case DISTANCE_Edit:
		{
		const float PctId = GetPctId(Seqi, Seqj, ColCount);
		if (PctId > 1.0)
			Quit("Internal error, DISTANCE_Edit, pct id=%.3g", PctId);
		d = (float) 1.0 - PctId;
		break;
		}
	default:
		Quit("DistCalcMSA: Invalid DISTANCE_%u", m_Distance);
		}
	return d;
	}

void DistCalcMSA::CalcDistRange(unsigned i, dist_t Dist[]) const
	{
	for (unsigned j = 0; j < i; ++j)
		Dist[j] = CalcDist(i, j);
	}

unsigned DistCalcMSA::GetCount() const
	{
	return m_ptrMSA->GetSeqCount();
	}

unsigned DistCalcMSA::GetId(unsigned i) const
	{
	return m_ptrMSA->GetSeqId(i);
	}

const char *DistCalcMSA::GetName(unsigned i) const
	{
	return m_ptrMSA->GetSeqName(i);
	}
