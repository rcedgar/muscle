#ifdef __MACH__

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/socket.h>
#include <sys/gmon.h>
#include <mach/vm_param.h>
#include <netinet/in.h>
#include <netinet/icmp6.h>
#include <sys/vmmeter.h>
#include <sys/proc.h>
#include <mach/task.h>
#include <mach/task_info.h>
#include <mach/mach_init.h>
#include <mach/vm_statistics.h>

const double DEFAULT_RAM = 1e9;
const double DEFAULT_MEM_USE = 1e6;

double GetNAN()
	{
	static unsigned long nan[2]={0xffffffff, 0x7fffffff};
	double dNAN = *( double* )nan;
	return dNAN;
	}

double g_dNAN = GetNAN();


double GetRAMSize()
	{
	static double CACHED_RAM = 0;
	if (CACHED_RAM != 0)
		return CACHED_RAM;

	uint64_t MemPages = 0;
	size_t Len = sizeof(MemPages);
	if (sysctlbyname("hw.memsize", &MemPages, &Len, NULL, 0) < 0)
		return DEFAULT_RAM;
	return (double) MemPages;
	}

double GetRAMSizeMB()
	{
	return GetRAMSize()/1e6;
	}

static double g_uPeakMemUseBytes;

double GetMaxMemUseBytes()
	{
	return g_uPeakMemUseBytes;
	}

double GetPeakMemUseBytes()
	{
	return GetMaxMemUseBytes();
	}

double GetMemUseBytes()
	{
	task_t mytask = mach_task_self();
	struct task_basic_info ti;
	memset((void *) &ti, 0, sizeof(ti));
	mach_msg_type_number_t count = TASK_BASIC_INFO_COUNT;
	kern_return_t ok = task_info(mytask, TASK_BASIC_INFO, (task_info_t) &ti, &count);
	if (ok == KERN_INVALID_ARGUMENT)
		return DEFAULT_MEM_USE;

	if (ok != KERN_SUCCESS)
		return DEFAULT_MEM_USE;

	double uBytes = (double ) ti.resident_size;
	if (uBytes > g_uPeakMemUseBytes)
		g_uPeakMemUseBytes = uBytes;
	return uBytes;
	}

double GetMemUseMB()
	{
	return GetMemUseBytes()/1e6;
	}

void OSInit()
	{
	}

#endif // __MACH__
