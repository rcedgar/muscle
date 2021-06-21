#include "muscle.h"
#include "edgelist.h"

EdgeList::EdgeList()
	{
	m_uNode1 = 0;
	m_uNode2 = 0;
	m_uCount = 0;
	m_uCacheSize = 0;
	}

EdgeList::~EdgeList()
	{
	Clear();
	}

void EdgeList::Clear()
	{
	delete[] m_uNode1;
	delete[] m_uNode2;
	m_uNode1 = 0;
	m_uNode2 = 0;
	m_uCount = 0;
	m_uCacheSize = 0;
	}

void EdgeList::Add(unsigned uNode1, unsigned uNode2)
	{
	if (m_uCount <= m_uCacheSize)
		Expand();
	m_uNode1[m_uCount] = uNode1;
	m_uNode2[m_uCount] = uNode2;
	++m_uCount;
	}

unsigned EdgeList::GetCount() const
	{
	return m_uCount;
	}

void EdgeList::GetEdge(unsigned uIndex, unsigned *ptruNode1, unsigned *ptruNode2) const
	{
	if (uIndex > m_uCount)
		Quit("EdgeList::GetEdge(%u) count=%u", uIndex, m_uCount);
	*ptruNode1 = m_uNode1[uIndex];
	*ptruNode2 = m_uNode2[uIndex];
	}

void EdgeList::Copy(const EdgeList &rhs)
	{
	Clear();
	const unsigned uCount = rhs.GetCount();
	for (unsigned n = 0; n < uCount; ++n)
		{
		unsigned uNode1;
		unsigned uNode2;
		rhs.GetEdge(n, &uNode1, &uNode2);
		Add(uNode1, uNode2);
		}
	}

void EdgeList::Expand()
	{
	unsigned uNewCacheSize = m_uCacheSize + 512;
	unsigned *NewNode1 = new unsigned[uNewCacheSize];
	unsigned *NewNode2 = new unsigned[uNewCacheSize];
	if (m_uCount > 0)
		{
		memcpy(NewNode1, m_uNode1, m_uCount*sizeof(unsigned));
		memcpy(NewNode2, m_uNode2, m_uCount*sizeof(unsigned));
		}
	delete[] m_uNode1;
	delete[] m_uNode2;
	m_uNode1 = NewNode1;
	m_uNode2 = NewNode2;
	m_uCacheSize = uNewCacheSize;
	}

void EdgeList::LogMe() const
	{
	for (unsigned n = 0; n < m_uCount; ++n)
		{
		if (n > 0)
			Log(" ");
		Log("%u->%u", m_uNode1[n], m_uNode2[n]);
		}
	Log("\n");
	}
