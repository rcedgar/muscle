#include "myutils.h"
#include "sort.h"
#include "sequence.h"
#include "multisequence.h"
#include "profpos3.h"
#include "profile3.h"
#include "objmgr.h"

static double g_SumAlloc;
static double g_LastSumAlloc;

void *myalloc_(size_t n, size_t m)
	{
	void *p = 0;
//#pragma omp critical
	{
	size_t Bytes = n*m;
	asserta(Bytes/n == m);

	static bool WarningDone;
	if (!WarningDone)
		{
		g_SumAlloc += Bytes;
		if (g_SumAlloc - g_LastSumAlloc > 1e9)
			{
			double RAM = GetPhysMemBytes();
			double Alloced = GetMemUseBytes();
			if (Alloced > RAM*0.9)
				{
				fprintf(stderr,
				  "\n\n"
				  "=========================================================\n"
				  "WARNING: %.3g Gb memory allocated so far.\n"
				  "This process may crash soon, or run slowly due to paging.\n"
				  "Typical cause of excessive memory use is large dataset\n"
				  "with long sequences (more than around 15k letters).\n"
				  "=========================================================\n"
				  "\n",
				  Alloced/1e9);
				WarningDone = true;
				}
			}
		g_LastSumAlloc = g_SumAlloc;
		}

	p = malloc(Bytes);
	if (p == 0)
		{
		double Alloced = GetMemUseBytes();
		Die("Out of memory mymalloc(%.3g), alloced %.3g bytes",
		  double(Bytes), double(Alloced));
		}
	}
	return p;
	}

void myfree_(void *p)
	{
	if (p == 0)
		return;
	free(p);
	}

#if TRACE_ALLOC
static map<void *, string> m_PtrToLocStr;
static map<void *, double> m_PtrToBytes;
static map<string, double> m_LocStrToBytes;

void *myalloc_track(const char *FileName, int LineNr, size_t n, size_t m)
	{
	void *p = myalloc_(n, m);
	double Bytes = double(n)*double(m);
	char Tmp[1024];
	const char *ptrName = FileName;
	for (const char *ptr = FileName; *ptr != 0; ++ptr)
		{
		char c = *ptr;
		if (c == '/' || c == '\\')
			ptrName = ptr + 1;
		}
	sprintf(Tmp, "%s:%d", ptrName, LineNr);
	string LocStr(Tmp);

#pragma omp critical
	{
	m_PtrToBytes[p] = Bytes;
	m_PtrToLocStr[p] = LocStr;
	if (m_LocStrToBytes.find(LocStr) == m_LocStrToBytes.end())
		m_LocStrToBytes[LocStr] = Bytes;
	else
		m_LocStrToBytes[LocStr] += Bytes;
	}

	return p;
	}

void myfree_track(void *p)
	{
#pragma omp critical
	{
	map<void *, string>::const_iterator iter1 = m_PtrToLocStr.find(p);
	map<void *, double>::const_iterator iter2 = m_PtrToBytes.find(p);
	double Bytes = 0;
	if (iter2 != m_PtrToBytes.end())
		Bytes = iter2->second;

	string LocStr;
	if (iter1 != m_PtrToLocStr.end())
		{
		LocStr = iter1->second;
		map<string, double>::iterator iter3 = m_LocStrToBytes.find(LocStr);
		if (iter3 != m_LocStrToBytes.end())
			iter3->second -= Bytes;
		}

	if (iter1 != m_PtrToLocStr.end())
		m_PtrToLocStr.erase(iter1);
	if (iter2 != m_PtrToBytes.end())
		m_PtrToBytes.erase(iter2);
	}

	myfree_(p);
	}

static double GetSizeFromStr(const string &s)
	{
	const char *ptr = strchr(s.c_str(), '=') + 1;
	double Size = atof(ptr);
	return Size;
	}

void LogAllocs()
	{
	vector<string> LocStrs;
	vector<double> Sizes;
	for (map<string, double>::const_iterator iter = m_LocStrToBytes.begin();
	  iter != m_LocStrToBytes.end(); ++iter)
		{
		const string &LocStr = iter->first;
		double Size = iter->second;
		Sizes.push_back(Size);
		LocStrs.push_back(LocStr);
		}
	const uint N = SIZE(LocStrs);
	uint *Order = (uint *) malloc(N*sizeof(uint));
	QuickSortOrderDesc(Sizes.data(), N, Order);
	uint MAXSHOW = 100;
	double SumLeak = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint k = Order[i];
		double Bytes = Sizes[k];
		if (Bytes == 0)
			continue;
		if (i < MAXSHOW)
			Log("%10.10s  %s\n", MemBytesToStr(Bytes), LocStrs[k].c_str());
		SumLeak += Bytes;
		}
	Sequence::LogNewDeleteCounts();
	Log("MultiSequence new %u, delete %u\n",
	  MultiSequence::m_NewCount, MultiSequence::m_DeleteCount);
	Log("Profile3 new %u, delete %u\n",
	  Profile3::m_NewCount, Profile3::m_DeleteCount);
	Log("ProfPos3 new %u, delete %u\n",
	  ProfPos3::m_NewCount, ProfPos3::m_DeleteCount);
	ProgressLog("\n****** TOTAL LEAK %s *******\n\n", 
	  MemBytesToStr(SumLeak));
	ObjMgr::LogGlobalStats();
	}
#else // TRACE_ALLOC
void LogAllocs() {}
#endif // TRACE_ALLOC
