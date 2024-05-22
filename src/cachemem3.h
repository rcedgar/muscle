#pragma once

#pragma once

class CacheMem3
	{
public:
	unsigned uCachePrefixCountB = 0;
	unsigned uCachePrefixCountA = 0;
	float *CacheMCurr = 0;
	float *CacheMNext = 0;
	float *CacheMPrev = 0;
	float *CacheDRow = 0;
	char **CacheTB = 0;

	CacheMem3()
		{
		uCachePrefixCountB = 0;
		uCachePrefixCountA = 0;
		CacheMCurr = 0;
		CacheMNext = 0;
		CacheMPrev = 0;
		CacheDRow = 0;
		CacheTB = 0;
		}

	~CacheMem3()
		{
		FreeMem();
		}

	void FreeMem()
		{
		myfree(CacheMCurr);
		myfree(CacheMNext);
		myfree(CacheMPrev);
		myfree(CacheDRow);
		for (unsigned i = 0; i < uCachePrefixCountA; ++i)
			myfree(CacheTB[i]);
		myfree(CacheTB);

		uCachePrefixCountB = 0;
		uCachePrefixCountA = 0;
		CacheMCurr = 0;
		CacheMNext = 0;
		CacheMPrev = 0;
		CacheDRow = 0;
		CacheTB = 0;
		}

public:
	void AllocCache(unsigned uPrefixCountA, unsigned uPrefixCountB)
		{
		if (uPrefixCountA <= uCachePrefixCountA && uPrefixCountB <= uCachePrefixCountB)
			return;

		FreeMem();

		uCachePrefixCountA = uPrefixCountA + 1024;
		uCachePrefixCountB = uPrefixCountB + 1024;

		CacheMCurr = myalloc(float, uCachePrefixCountB);
		CacheMNext = myalloc(float, uCachePrefixCountB);
		CacheMPrev = myalloc(float, uCachePrefixCountB);
		CacheDRow = myalloc(float, uCachePrefixCountB);

		CacheTB = myalloc(char *, uCachePrefixCountA);
		for (unsigned i = 0; i < uCachePrefixCountA; ++i)
			CacheTB[i] = myalloc(char, uCachePrefixCountB);
		}
	};
