#pragma once

#include <stdint.h>

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

typedef uint64_t TICKS;

static inline TICKS GetClockTicks()
	{
#ifdef _MSC_VER
	unsigned int aux;
	_mm_lfence();
	uint64_t t = __rdtscp(&aux);
	_mm_lfence();
	return t;
#else
	unsigned aux;
	_mm_lfence();
	uint64_t t = __rdtscp(&aux);
	_mm_lfence();
	return t;
#endif
	}
