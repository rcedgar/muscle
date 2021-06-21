#include "muscle.h"
#include "distfunc.h"
//#include <assert.h>

DistFunc::DistFunc()
	{
	m_Dists = 0;
	m_uCount = 0;
	m_uCacheCount = 0;
	m_Names = 0;
	m_Ids = 0;
	}

DistFunc::~DistFunc()
	{
	if (0 != m_Names)
		{
		for (unsigned i = 0; i < m_uCount; ++i)
			free(m_Names[i]);
		}
	delete[] m_Dists;
	delete[] m_Names;
	delete[] m_Ids;
	}

float DistFunc::GetDist(unsigned uIndex1, unsigned uIndex2) const
	{
	return m_Dists[VectorIndex(uIndex1, uIndex2)];
	}

unsigned DistFunc::GetCount() const
	{
	return m_uCount;
	}

void DistFunc::SetCount(unsigned uCount)
	{
	m_uCount = uCount;
	if (uCount <= m_uCacheCount)
		return;
	delete[] m_Dists;
	m_Dists = new float[VectorLength()];
	m_Names = new char *[m_uCount];
	m_Ids = new unsigned[m_uCount];
	m_uCacheCount = uCount;

	memset(m_Names, 0, m_uCount*sizeof(char *));
	memset(m_Ids, 0xff, m_uCount*sizeof(unsigned));
	memset(m_Dists, 0, VectorLength()*sizeof(float));
	}

void DistFunc::SetDist(unsigned uIndex1, unsigned uIndex2, float dDist)
	{
	m_Dists[VectorIndex(uIndex1, uIndex2)] = dDist;
	m_Dists[VectorIndex(uIndex2, uIndex1)] = dDist;
	}

unsigned DistFunc::VectorIndex(unsigned uIndex1, unsigned uIndex2) const
	{
	assert(uIndex1 < m_uCount && uIndex2 < m_uCount);
	return uIndex1*m_uCount + uIndex2;
	}

unsigned DistFunc::VectorLength() const
	{
	return m_uCount*m_uCount;
	}

void DistFunc::SetName(unsigned uIndex, const char szName[])
	{
	assert(uIndex < m_uCount);
	m_Names[uIndex] = strsave(szName);
	}

void DistFunc::SetId(unsigned uIndex, unsigned uId)
	{
	assert(uIndex < m_uCount);
	m_Ids[uIndex] = uId;
	}

const char *DistFunc::GetName(unsigned uIndex) const
	{
	assert(uIndex < m_uCount);
	return m_Names[uIndex];
	}

unsigned DistFunc::GetId(unsigned uIndex) const
	{
	assert(uIndex < m_uCount);
	return m_Ids[uIndex];
	}

void DistFunc::LogMe() const
	{
	Log("DistFunc::LogMe count=%u\n", m_uCount);
	Log("                     ");
	for (unsigned i = 0; i < m_uCount; ++i)
		Log(" %7u", i);
	Log("\n");

	Log("                     ");
	for (unsigned i = 0; i < m_uCount; ++i)
		Log(" %7.7s", m_Names[i] ? m_Names[i] : "");
	Log("\n");

	for (unsigned i = 0; i < m_uCount; ++i)
		{
		Log("%4u  %10.10s  :  ", i, m_Names[i] ? m_Names[i] : "");
		for (unsigned j = 0; j <= i; ++j)
			Log(" %7.4g", GetDist(i, j));
		Log("\n");
		}
	}
