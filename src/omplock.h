#ifndef omplock_h
#define omplock_h

static omp_lock_t g_Lock;

static bool omp_lock_init()
	{
	omp_init_lock(&g_Lock);
	return true;
	}
static bool omp_lock_init_done = omp_lock_init();

static inline void Lock()
	{
	omp_set_lock(&g_Lock);
	}

static inline void Unlock()
	{
	omp_unset_lock(&g_Lock);
	}

#define LOCK()		Lock()
#define UNLOCK()	Unlock()

#endif // omplock_h
