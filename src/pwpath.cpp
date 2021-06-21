#include "muscle.h"
#include "pwpath.h"
#include "seq.h"
#include "textfile.h"
#include "msa.h"

PWPath::PWPath()
	{
	m_uArraySize = 0;
	m_uEdgeCount = 0;
	m_Edges = 0;
	}

PWPath::~PWPath()
	{
	Clear();
	}

void PWPath::Clear()
	{
	delete[] m_Edges;
	m_Edges = 0;
	m_uArraySize = 0;
	m_uEdgeCount = 0;
	}

void PWPath::ExpandPath(unsigned uAdditionalEdgeCount)
	{
	PWEdge *OldPath = m_Edges;
	unsigned uEdgeCount = m_uArraySize + uAdditionalEdgeCount;

	m_Edges = new PWEdge[uEdgeCount];
	m_uArraySize = uEdgeCount;
	if (m_uEdgeCount > 0)
		memcpy(m_Edges, OldPath, m_uEdgeCount*sizeof(PWEdge));
	delete[] OldPath;
	}

void PWPath::AppendEdge(const PWEdge &Edge)
	{
	if (0 == m_uArraySize || m_uEdgeCount + 1 == m_uArraySize)
		ExpandPath(200);

	m_Edges[m_uEdgeCount] = Edge;
	++m_uEdgeCount;
	}

void PWPath::AppendEdge(char cType, unsigned uPrefixLengthA, unsigned uPrefixLengthB)
	{
	PWEdge e;
	e.uPrefixLengthA = uPrefixLengthA;
	e.uPrefixLengthB = uPrefixLengthB;
	e.cType = cType;
	AppendEdge(e);
	}

void PWPath::PrependEdge(const PWEdge &Edge)
	{
	if (0 == m_uArraySize || m_uEdgeCount + 1 == m_uArraySize)
		ExpandPath(1000);
	if (m_uEdgeCount > 0)
		memmove(m_Edges + 1, m_Edges, sizeof(PWEdge)*m_uEdgeCount);
	m_Edges[0] = Edge;
	++m_uEdgeCount;
	}

const PWEdge &PWPath::GetEdge(unsigned uEdgeIndex) const
	{
	assert(uEdgeIndex < m_uEdgeCount);
	return m_Edges[uEdgeIndex];
	}

void PWPath::Validate() const
	{
	const unsigned uEdgeCount = GetEdgeCount();
	if (0 == uEdgeCount)
		return;
	const PWEdge &FirstEdge = GetEdge(0);
	const PWEdge &LastEdge = GetEdge(uEdgeCount - 1);
	unsigned uStartA = FirstEdge.uPrefixLengthA;
	unsigned uStartB = FirstEdge.uPrefixLengthB;
	if (FirstEdge.cType != 'I')
		--uStartA;
	if (FirstEdge.cType != 'D')
		--uStartB;

	unsigned uPrefixLengthA = FirstEdge.uPrefixLengthA;
	unsigned uPrefixLengthB = FirstEdge.uPrefixLengthB;
	for (unsigned uEdgeIndex = 1; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = GetEdge(uEdgeIndex);
		switch (Edge.cType)
			{
		case 'M':
			if (uPrefixLengthA + 1 != Edge.uPrefixLengthA)
				Quit("PWPath::Validate MA %u", uPrefixLengthA);
			if (uPrefixLengthB + 1 != Edge.uPrefixLengthB)
				Quit("PWPath::Validate MB %u", uPrefixLengthB);
			++uPrefixLengthA;
			++uPrefixLengthB;
			break;
		case 'D':
			if (uPrefixLengthA + 1 != Edge.uPrefixLengthA)
				Quit("PWPath::Validate DA %u", uPrefixLengthA);
			if (uPrefixLengthB != Edge.uPrefixLengthB)
				Quit("PWPath::Validate DB %u", uPrefixLengthB);
			++uPrefixLengthA;
			break;
		case 'I':
			if (uPrefixLengthA != Edge.uPrefixLengthA)
				Quit("PWPath::Validate IA %u", uPrefixLengthA);
			if (uPrefixLengthB + 1 != Edge.uPrefixLengthB)
				Quit("PWPath::Validate IB %u", uPrefixLengthB);
			++uPrefixLengthB;
			break;
			}
		}
	}

void PWPath::LogMe() const
	{
	for (unsigned uEdgeIndex = 0; uEdgeIndex < GetEdgeCount(); ++uEdgeIndex)
		{
		const PWEdge &Edge = GetEdge(uEdgeIndex);
		if (uEdgeIndex > 0)
			Log(" ");
		Log("%c%d.%d",
		  Edge.cType,
		  Edge.uPrefixLengthA,
		  Edge.uPrefixLengthB);
		if ((uEdgeIndex > 0 && uEdgeIndex%10 == 0) ||
		 uEdgeIndex == GetEdgeCount() - 1)
			Log("\n");
		}
	}

void PWPath::Copy(const PWPath &Path)
	{
	Clear();
	const unsigned uEdgeCount = Path.GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		AppendEdge(Edge);
		}
	}

void PWPath::FromMSAPair(const MSA &msaA, const MSA &msaB)
	{
	const unsigned uColCount = msaA.GetColCount();
	if (uColCount != msaB.GetColCount())
		Quit("PWPath::FromMSAPair, lengths differ");

	Clear();

	unsigned uPrefixLengthA = 0;
	unsigned uPrefixLengthB = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bIsGapA = msaA.IsGapColumn(uColIndex);
		bool bIsGapB = msaB.IsGapColumn(uColIndex);

		PWEdge Edge;
		char cType;
		if (!bIsGapA && !bIsGapB)
			{
			cType = 'M';
			++uPrefixLengthA;
			++uPrefixLengthB;
			}
		else if (bIsGapA && !bIsGapB)
			{
			cType = 'I';
			++uPrefixLengthB;
			}
		else if (!bIsGapA && bIsGapB)
			{
			cType = 'D';
			++uPrefixLengthA;
			}
		else
			{
			assert(bIsGapB && bIsGapA);
			continue;
			}

		Edge.cType = cType;
		Edge.uPrefixLengthA = uPrefixLengthA;
		Edge.uPrefixLengthB = uPrefixLengthB;
		AppendEdge(Edge);
		}
	}

// Very similar to HMMPath::FromFile, should consolidate.
void PWPath::FromFile(TextFile &File)
	{
	Clear();
	char szToken[1024];
	File.GetTokenX(szToken, sizeof(szToken));
	if (0 != strcmp(szToken, "Path"))
		Quit("Invalid path file (Path)");

	File.GetTokenX(szToken, sizeof(szToken));
	if (0 != strcmp(szToken, "edges"))
		Quit("Invalid path file (edges)");

	File.GetTokenX(szToken, sizeof(szToken));
	if (!IsValidInteger(szToken))
		Quit("Invalid path file (edges value)");

	const unsigned uEdgeCount = (unsigned) atoi(szToken);
	unsigned uEdgeIndex = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
	// index
		File.GetTokenX(szToken, sizeof(szToken));
		if (!IsValidInteger(szToken))
			Quit("Invalid path file, invalid index '%s'", szToken);
		unsigned n = (unsigned) atoi(szToken);
		if (n != uEdgeIndex)
			Quit("Invalid path file, expecting edge %u got %u", uEdgeIndex, n);

	// type
		File.GetTokenX(szToken, sizeof(szToken));
		if (1 != strlen(szToken))
			Quit("Invalid path file, expecting state, got '%s'", szToken);
		const char cType = szToken[0];
		if ('M' != cType && 'D' != cType && cType != 'I' && 'S' != cType)
			Quit("Invalid path file, expecting state, got '%c'", cType);

	// prefix length A
		File.GetTokenX(szToken, sizeof(szToken));
		if (!IsValidInteger(szToken))
			Quit("Invalid path file, bad prefix length A '%s'", szToken);
		const unsigned uPrefixLengthA = (unsigned) atoi(szToken);

	// prefix length B
		File.GetTokenX(szToken, sizeof(szToken));
		if (!IsValidInteger(szToken))
			Quit("Invalid path file, bad prefix length B '%s'", szToken);
		const unsigned uPrefixLengthB = (unsigned) atoi(szToken);

		PWEdge Edge;
		Edge.cType = cType;
		Edge.uPrefixLengthA = uPrefixLengthA;
		Edge.uPrefixLengthB = uPrefixLengthB;
		AppendEdge(Edge);
		}
	File.GetTokenX(szToken, sizeof(szToken));
	if (0 != strcmp(szToken, "//"))
		Quit("Invalid path file (//)");
	}

void PWPath::ToFile(TextFile &File) const
	{
	const unsigned uEdgeCount = GetEdgeCount();

	File.PutString("Path\n");
	File.PutFormat("edges %u\n", uEdgeCount);
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = GetEdge(uEdgeIndex);
		File.PutFormat("%u %c %u %u\n",
		  uEdgeIndex,
		  Edge.cType,
		  Edge.uPrefixLengthA,
		  Edge.uPrefixLengthB);
		}
	File.PutString("//\n");
	}

void PWPath::AssertEqual(const PWPath &Path) const
	{
	const unsigned uEdgeCount = GetEdgeCount();
	if (uEdgeCount != Path.GetEdgeCount())
		{
		Log("PWPath::AssertEqual, this=\n");
		LogMe();
		Log("\nOther path=\n");
		Path.LogMe();
		Log("\n");
		Quit("PWPath::AssertEqual, Edge count different %u %u\n",
		  uEdgeCount, Path.GetEdgeCount());
		}

	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &e1 = GetEdge(uEdgeIndex);
		const PWEdge &e2 = Path.GetEdge(uEdgeIndex);
		if (e1.cType != e2.cType || e1.uPrefixLengthA != e2.uPrefixLengthA ||
		  e1.uPrefixLengthB != e2.uPrefixLengthB)
			{
			Log("PWPath::AssertEqual, this=\n");
			LogMe();
			Log("\nOther path=\n");
			Path.LogMe();
			Log("\n");
			Log("This edge %c%u.%u, other edge %c%u.%u\n",
			  e1.cType, e1.uPrefixLengthA, e1.uPrefixLengthB,
			  e2.cType, e2.uPrefixLengthA, e2.uPrefixLengthB);
			Quit("PWPath::AssertEqual, edge %u different\n", uEdgeIndex);
			}
		}
	}

bool PWPath::Equal(const PWPath &Path) const
	{
	const unsigned uEdgeCount = GetEdgeCount();
	if (uEdgeCount != Path.GetEdgeCount())
		return false;

	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &e1 = GetEdge(uEdgeIndex);
		const PWEdge &e2 = Path.GetEdge(uEdgeIndex);
		if (e1.cType != e2.cType || e1.uPrefixLengthA != e2.uPrefixLengthA ||
		  e1.uPrefixLengthB != e2.uPrefixLengthB)
			return false;
		}
	return true;
	}

unsigned PWPath::GetMatchCount() const
	{
	unsigned uMatchCount = 0;
	const unsigned uEdgeCount = GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &e = GetEdge(uEdgeIndex);
		if ('M' == e.cType)
			++uMatchCount;
		}
	return uMatchCount;
	}

unsigned PWPath::GetInsertCount() const
	{
	unsigned uInsertCount = 0;
	const unsigned uEdgeCount = GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &e = GetEdge(uEdgeIndex);
		if ('I' == e.cType)
			++uInsertCount;
		}
	return uInsertCount;
	}

unsigned PWPath::GetDeleteCount() const
	{
	unsigned uDeleteCount = 0;
	const unsigned uEdgeCount = GetEdgeCount();
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &e = GetEdge(uEdgeIndex);
		if ('D' == e.cType)
			++uDeleteCount;
		}
	return uDeleteCount;
	}

void PWPath::FromStr(const char Str[])
	{
	Clear();
	unsigned uPrefixLengthA = 0;
	unsigned uPrefixLengthB = 0;
	while (char c = *Str++)
		{
		switch (c)
			{
		case 'M':
			++uPrefixLengthA;
			++uPrefixLengthB;
			break;
		case 'D':
			++uPrefixLengthA;
			break;
		case 'I':
			++uPrefixLengthB;
			break;
		default:
			Quit("PWPath::FromStr, invalid state %c", c);
			}
		AppendEdge(c, uPrefixLengthA, uPrefixLengthB);
		}
	}
