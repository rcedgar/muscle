#include "muscle.h"

static const char usage_blob[] =
#include "help.h"
	;
const char *usage_txt[] = { usage_blob };
int g_n_usage_txt = 1;

void Usage(FILE *f)
	{
	PrintBanner(f);
	fputs(usage_blob, f);
	}
