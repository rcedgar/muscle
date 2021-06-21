#include "muscle.h"

#if	WIN32
#include <windows.h>
#include <crtdbg.h>
#include <psapi.h>
#include <float.h>
#include <stdio.h>

void DebugPrintf(const char *szFormat, ...)
	{
	va_list ArgList;
	char szStr[4096];

	va_start(ArgList, szFormat);
	vsprintf(szStr, szFormat, ArgList);

	OutputDebugString(szStr);
	}

double GetNAN()
	{
	static unsigned long nan[2]={0xffffffff, 0x7fffffff};
	double dNAN = *( double* )nan;
	assert(_isnan(dNAN));
	return dNAN;
	}

double g_dNAN = GetNAN();

void chkmem(const char szMsg[])
	{
	if (!_CrtCheckMemory())
		Quit("chkmem(%s)", szMsg);
	}

void Break()
	{
	if (IsDebuggerPresent())
		DebugBreak();
	}

const char *GetCmdLine()
	{
	return GetCommandLine();
	}

static unsigned uPeakMemUseBytes;

double GetRAMSizeMB()
	{
	MEMORYSTATUS MS;
	GlobalMemoryStatus(&MS);
	return MS.dwAvailPhys/1e6;
	}

double GetMemUseMB()
	{
	HANDLE hProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS PMC;
	BOOL bOk = GetProcessMemoryInfo(hProc, &PMC, sizeof(PMC));
	assert(bOk);
	//printf("GetMemUseMB()\n");
	//printf("%12u  PageFaultCount\n", (unsigned) PMC.PageFaultCount);
	//printf("%12u  PagefileUsage\n", (unsigned) PMC.PagefileUsage);
	//printf("%12u  PeakPagefileUsage\n", (unsigned) PMC.PeakPagefileUsage);
	//printf("%12u  WorkingSetSize\n", (unsigned) PMC.WorkingSetSize);
	//printf("%12u  PeakWorkingSetSize\n", (unsigned) PMC.PeakWorkingSetSize);
	//printf("%12u  QuotaPagedPoolUsage\n", (unsigned) PMC.QuotaPagedPoolUsage);
	//printf("%12u  QuotaPeakPagedPoolUsage\n", (unsigned) PMC.QuotaPeakPagedPoolUsage);
	//printf("%12u  QuotaNonPagedPoolUsage\n", (unsigned) PMC.QuotaNonPagedPoolUsage);
	//printf("%12u  QuotaPeakNonPagedPoolUsage\n", (unsigned) PMC.QuotaPeakNonPagedPoolUsage);
	unsigned uBytes = (unsigned) PMC.WorkingSetSize;
	if (uBytes > uPeakMemUseBytes)
		uPeakMemUseBytes = uBytes;
	return (uBytes + 500000.0)/1000000.0;
	}

double GetPeakMemUseMB()
	{
	return (uPeakMemUseBytes + 500000.0)/1000000.0;
	}

void CheckMemUse()
	{
// Side-effect: sets peak usage in uPeakMemUseBytes
	GetMemUseMB();
	}

double GetCPUGHz()
	{
	double dGHz = 2.5;
	const char *e = getenv("CPUGHZ");
	if (0 != e)
		dGHz = atof(e);
	if (dGHz < 0.1 || dGHz > 1000.0)
		Quit("Invalid value '%s' for environment variable CPUGHZ", e);
	return dGHz;
	}
#endif	// WIN32
