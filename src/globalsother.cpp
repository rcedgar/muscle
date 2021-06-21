#include "muscle.h"

#if		!defined(__linux__) && !defined(_MSC_VER) && !defined(__MACH__)

double GetNAN()
	{
	return 0.0;
	}

double g_dNAN = GetNAN();

void chkmem(const char szMsg[])
	{
	}

void Break()
	{
	}

char szCmdLine[4096];

const char *GetCmdLine()
	{
	return "muscle";
	}

double GetMemUseMB()
	{
	return 100.0;
	}

void SaveCmdLine(int argc, char *argv[])
	{
	for (int i = 0; i < argc; ++i)
		{
		if (i > 0)
			strcat(szCmdLine, " ");
		strcat(szCmdLine, argv[i]);
		}
	}

double GetPeakMemUseMB()
	{
	return 100.0;
	}

double GetCPUGHz()
	{
	return 2.0;
	}

void CheckMemUse()
	{
	}

double GetRAMSizeMB()
	{
	return 500.0;
	}

#endif

