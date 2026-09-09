#include "myutils.h"
#include "git_hash.h"

#define FLAG_OPT(Name)	bool opt_##Name; bool optset_##Name; bool optused_##Name;
#define UNS_OPT(Name)	unsigned opt_##Name; bool optset_##Name; bool optused_##Name;
#define FLT_OPT(Name)	double opt_##Name; bool optset_##Name; bool optused_##Name;
#define STR_OPT(Name)	const char *opt_##Name = ""; bool optset_##Name; bool optused_##Name;
#include "myopts.h"

const char *g_ProgramName = PROGRAM_NAME;
const char *g_MyVersion = MY_VERSION;
const char *g_GitHash = GIT_HASH;
