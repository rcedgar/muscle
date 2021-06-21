#ifndef enumopts_h
#define enumopts_h

struct EnumOpt
	{
	const char *pstrOpt;
	int iValue;
	};

#define	s(t)		extern EnumOpt t##_Opts[];
#define c(t, x)		/* empty */
#define e(t)		/* empty */
#include "enums.h"	


#endif // enumopts_h
