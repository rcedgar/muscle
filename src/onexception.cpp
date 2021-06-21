#include "muscle.h"
#include <stdio.h>

static char szOnExceptionMessage[] =
	{
	"\nFatal error, exception caught.\n"
	};

void OnException()
	{
	fprintf(stderr, "%s", szOnExceptionMessage);
	Log("%s", szOnExceptionMessage);
	Log("Finished %s\n", GetTimeAsStr());
	exit(EXIT_Except);
	}
