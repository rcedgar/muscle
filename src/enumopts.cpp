#include "muscle.h"
#include "enumopts.h"

#define	s(t)		EnumOpt t##_Opts[] = {
#define c(t, x)		#x, t##_##x,
#define e(t)		0, 0 };

#include "enums.h"
