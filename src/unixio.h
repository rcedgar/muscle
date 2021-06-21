#ifdef	WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

#if	!defined(WIN32) && !defined(O_BINARY)
#define	O_BINARY	0
#endif
