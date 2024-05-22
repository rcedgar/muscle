#ifndef objtype_h
#define objtype_h

enum ObjType
	{
#define T(x)	OT_##x,
#include "objtypes.h"
	OTCount
	};

#endif // objtype_h
