#include "muscle.h"

#if		defined(__linux__)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>

const int ONE_MB = 1000000;
const int MEM_WARNING_THRESHOLD = 20*ONE_MB;

double GetNAN()
	{
	static unsigned long nan[2]={0xffffffff, 0x7fffffff};
	double dNAN = *( double* )nan;
	return dNAN;
	}

double g_dNAN = GetNAN();

void chkmem(const char szMsg[])
	{
	//assert(_CrtCheckMemory());
	}

void Break()
	{
	//DebugBreak();
	}

static char szCmdLine[4096];

void *ptrStartBreak = sbrk(0);

const char *GetCmdLine()
	{
	return szCmdLine;
	}

double GetMemUseMB()
	{
	static char statm[64];
	static int PageSize;
	if (0 == statm[0])
		{
		PageSize = sysconf(_SC_PAGESIZE);
		pid_t pid = getpid();
		sprintf(statm, "/proc/%d/statm", (int) pid);
		}

	int fd = open(statm, O_RDONLY);
	if (-1 == fd)
		return -1;
	char Buffer[64];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		{
		static bool Warned = false;
		if (!Warned)
			{
			Warned = true;
			Warning("*Warning* Cannot read %s errno=%d %s",
			  statm, errno, strerror(errno));
			}
		return 0;
		}
	Buffer[n] = 0;
	int Pages = atoi(Buffer);

	return ((double) Pages * (double) PageSize)/1e6;
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

double dPeakMemUseMB = 0;

double GetPeakMemUseMB()
	{
	CheckMemUse();
	return dPeakMemUseMB;
	}

double GetCPUGHz()
	{
	double dGHz = 2.5;
	const char *e = getenv("CPUGHZ");
	if (0 != e)
		dGHz = atof(e);
	return dGHz;
	}

void CheckMemUse()
	{
	double dMB = GetMemUseMB();
	if (dMB > dPeakMemUseMB)
		dPeakMemUseMB = dMB;
	}

double GetRAMSizeMB()
	{
	const double DEFAULT_RAM = 500;
	static double RAMMB = 0;
	if (RAMMB != 0)
		return RAMMB;

	int fd = open("/proc/meminfo", O_RDONLY);
	if (-1 == fd)
		{
		static bool Warned = false;
		if (!Warned)
			{
			Warned = true;
			Warning("*Warning* Cannot open /proc/meminfo errno=%d %s",
			  errno, strerror(errno));
			}
		return DEFAULT_RAM;
		}
	char Buffer[1024];
	int n = read(fd, Buffer, sizeof(Buffer) - 1);
	close(fd);
	fd = -1;

	if (n <= 0)
		{
		static bool Warned = false;
		if (!Warned)
			{
			Warned = true;
			Warning("*Warning* Cannot read /proc/meminfo errno=%d %s",
			  errno, strerror(errno));
			}
		return DEFAULT_RAM;
		}
	Buffer[n] = 0;
	char *pMem = strstr(Buffer, "MemTotal: ");
	if (0 == pMem)
		{
		static bool Warned = false;
		if (!Warned)
			{
			Warned = true;
			Warning("*Warning* 'MemTotal:' not found in /proc/meminfo");
			}
		return DEFAULT_RAM;
		}
	int Bytes = atoi(pMem+9)*1000;
	return ((double) Bytes)/1e6;
	}

#endif	// !WIN32
