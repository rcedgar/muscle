#include "muscle.h"
#include "pwpath.h"
#include "estring.h"
#include "seq.h"
#include "msa.h"

/***
An "estring" is an edit string that operates on a sequence.
An estring is represented as a vector of integers.
It is interpreted in order of increasing suffix.
A positive value n means copy n letters.
A negative value -n means insert n indels.
Zero marks the end of the vector.
Consecutive entries must have opposite sign, i.e. the
shortest possible representation must be used.

A "tpair" is a traceback path for a pairwise alignment
represented as two estrings, one for each sequence.
***/

#define c2(c,d)	(((unsigned char) c) << 8 | (unsigned char) d)

unsigned LengthEstring(const int es[])
	{
	unsigned i = 0;
	while (*es++ != 0)
		++i;
	return i;
	}

int *EstringNewCopy(const int es[])
	{
	unsigned n = LengthEstring(es) + 1;
	int *esNew = new int[n];
	memcpy(esNew, es, n*sizeof(int));
	return esNew;
	}

void LogEstring(const int es[])
	{
	Log("<");
	for (unsigned i = 0; es[i] != 0; ++i)
		{
		if (i > 0)
			Log(" ");
		Log("%d", es[i]);
		}
	Log(">");
	}

static bool EstringsEq(const int es1[], const int es2[])
	{
	for (;;)
		{
		if (*es1 != *es2)
			return false;
		if (0 == *es1)
			break;
		++es1;
		++es2;
		}
	return true;
	}

static void EstringCounts(const int es[], unsigned *ptruSymbols,
  unsigned *ptruIndels)
	{
	unsigned uSymbols = 0;
	unsigned uIndels = 0;
	for (unsigned i = 0; es[i] != 0; ++i)
		{
		int n = es[i];
		if (n > 0)
			uSymbols += n;
		else if (n < 0)
			uIndels += -n;
		}
	*ptruSymbols = uSymbols;
	*ptruIndels = uIndels;
	}

static char *EstringOp(const int es[], const char s[])
	{
	unsigned uSymbols;
	unsigned uIndels;
	EstringCounts(es, &uSymbols, &uIndels);
	assert((unsigned) strlen(s) == uSymbols);
	char *sout = new char[uSymbols + uIndels + 1];
	char *psout = sout;
	for (;;)
		{
		int n = *es++;
		if (0 == n)
			break;
		if (n > 0)
			for (int i = 0; i < n; ++i)
				*psout++ = *s++;
		else
			for (int i = 0; i < -n; ++i)
				*psout++ = '-';
		}
	assert(0 == *s);
	*psout = 0;
	return sout;
	}

void EstringOp(const int es[], const Seq &sIn, Seq &sOut)
	{
#if	DEBUG
	unsigned uSymbols;
	unsigned uIndels;
	EstringCounts(es, &uSymbols, &uIndels);
	assert(sIn.Length() == uSymbols);
#endif
	sOut.Clear();
	sOut.SetName(sIn.GetName());
	int p = 0;
	for (;;)
		{
		int n = *es++;
		if (0 == n)
			break;
		if (n > 0)
			for (int i = 0; i < n; ++i)
				{
				char c = sIn[p++];
				sOut.push_back(c);
				}
		else
			for (int i = 0; i < -n; ++i)
				sOut.push_back('-');
		}
	}

unsigned EstringOp(const int es[], const Seq &sIn, MSA &a)
	{
	unsigned uSymbols;
	unsigned uIndels;
	EstringCounts(es, &uSymbols, &uIndels);
	assert(sIn.Length() == uSymbols);

	unsigned uColCount = uSymbols + uIndels;

	a.Clear();
	a.SetSize(1, uColCount);

	a.SetSeqName(0, sIn.GetName());
	a.SetSeqId(0, sIn.GetId());

	unsigned p = 0;
	unsigned uColIndex = 0;
	for (;;)
		{
		int n = *es++;
		if (0 == n)
			break;
		if (n > 0)
			for (int i = 0; i < n; ++i)
				{
				char c = sIn[p++];
				a.SetChar(0, uColIndex++, c);
				}
		else
			for (int i = 0; i < -n; ++i)
				a.SetChar(0, uColIndex++, '-');
		}
	assert(uColIndex == uColCount);
	return uColCount;
	}

void PathToEstrings(const PWPath &Path, int **ptresA, int **ptresB)
	{
// First pass to determine size of estrings esA and esB
	const unsigned uEdgeCount = Path.GetEdgeCount();
	if (0 == uEdgeCount)
		{
		int *esA = new int[1];
		int *esB = new int[1];
		esA[0] = 0;
		esB[0] = 0;
		*ptresA = esA;
		*ptresB = esB;
		return;
		}

	unsigned iLengthA = 1;
	unsigned iLengthB = 1;
	const char cFirstEdgeType = Path.GetEdge(0).cType;
	char cPrevEdgeType = cFirstEdgeType;
	for (unsigned uEdgeIndex = 1; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		char cEdgeType = Edge.cType;

		switch (c2(cPrevEdgeType, cEdgeType))
			{
		case c2('M', 'M'):
		case c2('D', 'D'):
		case c2('I', 'I'):
			break;

		case c2('D', 'M'):
		case c2('M', 'D'):
			++iLengthB;
			break;

		case c2('I', 'M'):
		case c2('M', 'I'):
			++iLengthA;
			break;

		case c2('I', 'D'):
		case c2('D', 'I'):
			++iLengthB;
			++iLengthA;
			break;

		default:
			assert(false);
			}
		cPrevEdgeType = cEdgeType;
		}

// Pass2 for seq A
	{
	int *esA = new int[iLengthA+1];
	unsigned iA = 0;
	switch (Path.GetEdge(0).cType)
		{
	case 'M':
	case 'D':
		esA[0] = 1;
		break;

	case 'I':
		esA[0] = -1;
		break;

	default:
		assert(false);
		}

	char cPrevEdgeType = cFirstEdgeType;
	for (unsigned uEdgeIndex = 1; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		char cEdgeType = Edge.cType;

		switch (c2(cPrevEdgeType, cEdgeType))
			{
		case c2('M', 'M'):
		case c2('D', 'D'):
		case c2('D', 'M'):
		case c2('M', 'D'):
			++(esA[iA]);
			break;

		case c2('I', 'D'):
		case c2('I', 'M'):
			++iA;
			esA[iA] = 1;
			break;

		case c2('M', 'I'):
		case c2('D', 'I'):
			++iA;
			esA[iA] = -1;
			break;

		case c2('I', 'I'):
			--(esA[iA]);
			break;

		default:
			assert(false);
			}

		cPrevEdgeType = cEdgeType;
		}
	assert(iA == iLengthA - 1);
	esA[iLengthA] = 0;
	*ptresA = esA;
	}

	{
// Pass2 for seq B
	int *esB = new int[iLengthB+1];
	unsigned iB = 0;
	switch (Path.GetEdge(0).cType)
		{
	case 'M':
	case 'I':
		esB[0] = 1;
		break;

	case 'D':
		esB[0] = -1;
		break;

	default:
		assert(false);
		}

	char cPrevEdgeType = cFirstEdgeType;
	for (unsigned uEdgeIndex = 1; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		char cEdgeType = Edge.cType;

		switch (c2(cPrevEdgeType, cEdgeType))
			{
		case c2('M', 'M'):
		case c2('I', 'I'):
		case c2('I', 'M'):
		case c2('M', 'I'):
			++(esB[iB]);
			break;

		case c2('D', 'I'):
		case c2('D', 'M'):
			++iB;
			esB[iB] = 1;
			break;

		case c2('M', 'D'):
		case c2('I', 'D'):
			++iB;
			esB[iB] = -1;
			break;

		case c2('D', 'D'):
			--(esB[iB]);
			break;

		default:
			assert(false);
			}

		cPrevEdgeType = cEdgeType;
		}
	assert(iB == iLengthB - 1);
	esB[iLengthB] = 0;
	*ptresB = esB;
	}

#if	DEBUG
	{
	const PWEdge &LastEdge = Path.GetEdge(uEdgeCount - 1);
	unsigned uSymbols;
	unsigned uIndels;
	EstringCounts(*ptresA, &uSymbols, &uIndels);
	assert(uSymbols == LastEdge.uPrefixLengthA);
	assert(uSymbols + uIndels == uEdgeCount);

	EstringCounts(*ptresB, &uSymbols, &uIndels);
	assert(uSymbols == LastEdge.uPrefixLengthB);
	assert(uSymbols + uIndels == uEdgeCount);

	PWPath TmpPath;
	EstringsToPath(*ptresA, *ptresB, TmpPath);
	TmpPath.AssertEqual(Path);
	}
#endif
	}

void EstringsToPath(const int esA[], const int esB[], PWPath &Path)
	{
	Path.Clear();
	unsigned iA = 0;
	unsigned iB = 0;
	int nA = esA[iA++];
	int nB = esB[iB++];
	unsigned uPrefixLengthA = 0;
	unsigned uPrefixLengthB = 0;
	for (;;)
		{
		char cType;
		if (nA > 0)
			{
			if (nB > 0)
				{
				cType = 'M';
				--nA;
				--nB;
				}
			else if (nB < 0)
				{
				cType = 'D';
				--nA;
				++nB;
				}
			else
				assert(false);
			}
		else if (nA < 0)
			{
			if (nB > 0)
				{
				cType = 'I';
				++nA;
				--nB;
				}
			else
				assert(false);
			}
		else
			assert(false);

		switch (cType)
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
			}

		PWEdge Edge;
		Edge.cType = cType;
		Edge.uPrefixLengthA = uPrefixLengthA;
		Edge.uPrefixLengthB = uPrefixLengthB;
		Path.AppendEdge(Edge);

		if (nA == 0)
			{
			if (0 == esA[iA])
				{
				assert(0 == esB[iB]);
				break;
				}
			nA = esA[iA++];
			}
		if (nB == 0)
			nB = esB[iB++];
		}
	}

/***
Multiply two estrings to make a third estring.
The product of two estrings e1*e2 is defined to be
the estring that produces the same result as applying
e1 then e2. Multiplication is not commutative. In fact,
the reversed order is undefined unless both estrings
consist of a single, identical, positive entry.
A primary motivation for using estrings is that
multiplication is very fast, reducing the time
needed to construct the root alignment.

Example

	<-1,3>(XXX)	= -XXX
	<2,-1,2>(-XXX) = -X-XX

Therefore,

	<-1,3>*<2,-1,2> = <-1,1,-1,2>
***/

static bool CanMultiplyEstrings(const int es1[], const int es2[])
	{
	unsigned uSymbols1;
	unsigned uSymbols2;
	unsigned uIndels1;
	unsigned uIndels2;
	EstringCounts(es1, &uSymbols1, &uIndels1);
	EstringCounts(es2, &uSymbols2, &uIndels2);
	return uSymbols1 + uIndels1 == uSymbols2;
	}

static inline void AppendGaps(int esp[], int &ip, int n)
	{
	assert(n < SHRT_MAX);
	if (-1 == ip)
		esp[++ip] = n;
	else if (esp[ip] < 0)
		esp[ip] += n;
	else
		esp[++ip] = n;
	}

static inline void AppendSymbols(int esp[], int &ip, int n)
	{
	assert(n < SHRT_MAX);
	if (-1 == ip)
		esp[++ip] = n;
	else if (esp[ip] > 0)
		esp[ip] += n;
	else
		esp[++ip] = n;
	}

void MulEstrings(const int es1[], const int es2[], int esp[])
	{
	unsigned i1 = 0;
	int ip = -1;
	int n1 = es1[i1++];
	for (unsigned i2 = 0; ; ++i2)
		{
		int n2 = es2[i2];
		if (0 == n2)
			break;
		if (n2 > 0)
			{
			for (;;)
				{
				if (n1 < 0)
					{
					if (n2 > -n1)
						{
						AppendGaps(esp, ip, n1);
						n2 += n1;
						n1 = es1[i1++];
						}
					else if (n2 == -n1)
						{
						AppendGaps(esp, ip, n1);
						n1 = es1[i1++];
						break;
						}
					else
						{
						assert(n2 < -n1);
						AppendGaps(esp, ip, -n2);
						n1 += n2;
						break;
						}
					}
				else
					{
					assert(n1 > 0);
					if (n2 > n1)
						{
						AppendSymbols(esp, ip, n1);
						n2 -= n1;
						n1 = es1[i1++];
						}
					else if (n2 == n1)
						{
						AppendSymbols(esp, ip, n1);
						n1 = es1[i1++];
						break;
						}
					else
						{
						assert(n2 < n1);
						AppendSymbols(esp, ip, n2);
						n1 -= n2;
						break;
						}
					}
				}
			}
		else
			{
			assert(n2 < 0);
			AppendGaps(esp, ip, n2);
			}
		}
	esp[++ip] = 0;

#if	DEBUG
	{
	int MaxLen = (int) (LengthEstring(es1) + LengthEstring(es2) + 1);
	assert(ip < MaxLen);
	if (ip >= 2)
		for (int i = 0; i < ip - 2; ++i)
			{
			if (!(esp[i] > 0 && esp[i+1] < 0 || esp[i] < 0 && esp[i+1] > 0))
				{
				Log("Bad result of MulEstring: ");
				LogEstring(esp);
				Quit("Assert failed (alternating signs)");
				}
			}
	unsigned uSymbols1;
	unsigned uSymbols2;
	unsigned uSymbolsp;
	unsigned uIndels1;
	unsigned uIndels2;
	unsigned uIndelsp;
	EstringCounts(es1, &uSymbols1, &uIndels1);
	EstringCounts(es2, &uSymbols2, &uIndels2);
	EstringCounts(esp, &uSymbolsp, &uIndelsp);
	if (uSymbols1 + uIndels1 != uSymbols2)
		{
		Log("Bad result of MulEstring: ");
		LogEstring(esp);
		Quit("Assert failed (counts1 %u %u %u)",
		  uSymbols1, uIndels1, uSymbols2);
		}
	}
#endif
	}

static void test(const int es1[], const int es2[], const int esa[])
	{
	unsigned uSymbols1;
	unsigned uSymbols2;
	unsigned uIndels1;
	unsigned uIndels2;
	EstringCounts(es1, &uSymbols1, &uIndels1);
	EstringCounts(es2, &uSymbols2, &uIndels2);

	char s[4096];
	memset(s, 'X', sizeof(s));
	s[uSymbols1] = 0;

	char *s1 = EstringOp(es1, s);
	char *s12 = EstringOp(es2, s1);

	memset(s, 'X', sizeof(s));
	s[uSymbols2] = 0;
	char *s2 = EstringOp(es2, s);

	Log("%s * %s = %s\n", s1, s2, s12);

	LogEstring(es1);
	Log(" * ");
	LogEstring(es2);
	Log(" = ");
	LogEstring(esa);
	Log("\n");

	int esp[4096];
	MulEstrings(es1, es2, esp);
	LogEstring(esp);
	if (!EstringsEq(esp, esa))
		Log(" *ERROR* ");
	Log("\n");

	memset(s, 'X', sizeof(s));
	s[uSymbols1] = 0;
	char *sp = EstringOp(esp, s);
	Log("%s\n", sp);
	Log("\n==========\n\n");
	}

void TestEstrings()
	{
	SetListFileName("c:\\tmp\\muscle.log", false);
	//{
	//int es1[] = { -1, 1, -1, 0 };
	//int es2[] = { 1, -1, 2, 0 };
	//int esa[] = { -2, 1, -1, 0 };
	//test(es1, es2, esa);
	//}
	//{
	//int es1[] = { 2, -1, 2, 0 };
	//int es2[] = { 1, -1, 3, -1, 1, 0 };
	//int esa[] = { 1, -1, 1, -1, 1, -1, 1, 0 };
	//test(es1, es2, esa);
	//}
	//{
	//int es1[] = { -1, 3, 0 };
	//int es2[] = { 2, -1, 2, 0 };
	//int esa[] = { -1, 1, -1, 2, 0 };
	//test(es1, es2, esa);
	//}
	//{
	//int es1[] = { -1, 1, -1, 1, 0};
	//int es2[] = { 4, 0 };
	//int esa[] = { -1, 1, -1, 1, 0};
	//test(es1, es2, esa);
	//}
	//{
	//int es1[] = { 1, -1, 1, -1, 0};
	//int es2[] = { 4, 0 };
	//int esa[] = { 1, -1, 1, -1, 0};
	//test(es1, es2, esa);
	//}
	//{
	//int es1[] = { 1, -1, 1, -1, 0};
	//int es2[] = { -1, 4, -1, 0 };
	//int esa[] = { -1, 1, -1, 1, -2, 0};
	//test(es1, es2, esa);
	//}
	{
	int es1[] = { 106, -77, 56, -2, 155, -3, 123, -2, 0};
	int es2[] = { 50, -36, 34, -3, 12, -6, 1, -6, 18, -17, 60, -5, 349, -56, 0 };
	int esa[] = { 0 };
	test(es1, es2, esa);
	}
	exit(0);
	}
