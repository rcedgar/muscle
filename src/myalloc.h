#pragma once

#define TRACE_ALLOC	0

void *myalloc_(size_t n, size_t m);
void myfree_(void *p);
void *myalloc_track(const char *FileName, int LineNr, size_t n, size_t m);
void myfree_track(void *p);
void LogAllocs();

#if TRACE_ALLOC
#define myalloc(t, m)	(t *) myalloc_track(__FILE__, __LINE__, sizeof(t), (m))
#define myfree(p)		myfree_track(p)
#else
#define myalloc(t, m)	(t *) myalloc_(sizeof(t), (m))
#define myfree(p)		myfree_(p)
#endif
