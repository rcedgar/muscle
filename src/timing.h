#ifndef getticks_h
#define getticks_h

#if 0

// ~3 x 10^9 ticks/sec

#ifdef _MSC_VER
#include <Windows.h>
typedef unsigned __int64 TICKS;

#define	GetClockTicks	__rdtsc

#elif defined(__APPLE__)
typedef uint64_t TICKS;
__inline__ uint64_t GetClockTicks()
	{
	return 0;
	}

#elif __GNUC__

typedef uint64_t TICKS;
__inline__ uint64_t GetClockTicks()
	{
	uint32_t lo, hi;
	/* We cannot use "=A", since this would use %rax on x86_64 */
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (uint64_t)hi << 32 | lo;
	}

#else
#error	"getticks_h, unknown compiler"
#endif

#endif

#endif // getticks_h
