#include "muscle.h"

void Usage(FILE *f)
	{
	PrintBanner(f);
	fputs(
#include "help.h"
	, f);
	}
