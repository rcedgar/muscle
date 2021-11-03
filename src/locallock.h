#pragma once

#include <omp.h>

static omp_lock_t g_Lock;
static bool InitLock()
	{
	omp_init_lock(&g_Lock);
	return true;
	}
static bool g_InitDone = InitLock();

static void Lock()
	{
	omp_set_lock(&g_Lock);
	}

static void Unlock()
	{
	omp_unset_lock(&g_Lock);
	}
