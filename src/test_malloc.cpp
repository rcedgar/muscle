#include "myutils.h"

static void Test(size_t n)
	{
	void *p = malloc(n);
	const char *s = MemBytesToStr(n);
	ProgressLog("%10.10s  %s\n", s, (p == 0 ? "failed" : "ok"));
	if (p != 0)
		free(p);
	}

void cmd_test_malloc()
	{
	opt(test_malloc);

	Test(18398178000);

	size_t n = (uint64(UINT_MAX) + 1)/2;
	for (uint i = 0; i < 8; ++i)
		{
		Test(n);
		n *= 2;
		}
	}
