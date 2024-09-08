#include "myutils.h"
#include "enumpaths.h"

uint GetNA(const string &Path)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = Path[i];
		if (c == 'M' || c == 'D')
			++n;
		}
	return n;
	}

uint GetNB(const string &Path)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = Path[i];
		if (c == 'M' || c == 'I')
			++n;
		}
	return n;
	}

void EnumPathsLocalRecurse(uint LoA, uint HiA, uint LoB, uint HiB,
  string Path, OnPathLocal_fn OnPath)
	{
	const uint ColCount = SIZE(Path);
	asserta(ColCount > 0);
	asserta(Path[0] == 'M');
	asserta(LoA <= HiA);
	asserta(LoB <= HiB);
	uint SubLA = HiA - LoA + 1;
	uint SubLB = HiB - LoB + 1;
	uint NA = GetNA(Path);
	uint NB = GetNB(Path);
	asserta(NA <= SubLA && NB <= SubLB);
	if (NA == SubLA && NB == SubLB && Path[ColCount-1] == 'M')
		{
		OnPath(LoA, LoB, Path);
		return;
		}
	char c = Path[ColCount-1];
	if (NA < SubLA && NB < SubLB)
		EnumPathsLocalRecurse(LoA, HiA, LoB, HiB, Path + 'M', OnPath);

	if (NA < SubLA && NB <= SubLB && c != 'I')
		EnumPathsLocalRecurse(LoA, HiA, LoB, HiB, Path + 'D', OnPath);

	if (NA <= SubLA && NB < SubLB && c != 'D')
		EnumPathsLocalRecurse(LoA, HiA, LoB, HiB, Path + 'I', OnPath);
	}

void EnumPathsLocal(uint LoA, uint HiA, uint LA, 
  uint LoB, uint HiB, uint LB, OnPathLocal_fn OnPath)
	{
	EnumPathsLocalRecurse(LoA, HiA, LoB, HiB, "M", OnPath);
	}

void EnumPathsLocal(uint LA, uint LB, OnPathLocal_fn OnPath)
	{
	for (uint LoA = 0; LoA < LA; ++LoA)
		for (uint HiA = LoA; HiA < LA; ++HiA)
			for (uint LoB = 0; LoB < LB; ++LoB)
				for (uint HiB = LoB; HiB < LB; ++HiB)
					EnumPathsLocal(LoA, HiA, LA, LoB, HiB, LB, OnPath);
	}

void EnumPathsGlobalRecurse(uint LA, uint LB,
  string Path, OnPathGlobal_fn OnPath)
	{
	const uint ColCount = SIZE(Path);
	asserta(ColCount > 0);
	uint NA = GetNA(Path);
	uint NB = GetNB(Path);
	asserta(NA <= LA && NB <= LB);
	if (NA == LA && NB == LB)
		{
		OnPath(Path);
		return;
		}
	char c = Path[ColCount-1];
	if (NA < LA && NB < LB)
		EnumPathsGlobalRecurse(LA, LB, Path + 'M', OnPath);

	if (NA < LA && NB <= LB && c != 'I')
		EnumPathsGlobalRecurse(LA, LB, Path + 'D', OnPath);

	if (NA <= LA && NB < LB && c != 'D')
		EnumPathsGlobalRecurse(LA, LB, Path + 'I', OnPath);
	}

void EnumPathsGlobal(uint LA, uint LB, OnPathGlobal_fn OnPath)
	{
	EnumPathsGlobalRecurse(LA, LB, "M", OnPath);
	EnumPathsGlobalRecurse(LA, LB, "D", OnPath);
	EnumPathsGlobalRecurse(LA, LB, "I", OnPath);
	}

void OnPathGlobal(const string &Path)
	{
	Log("%s\n", Path.c_str());
	}

void OnPathLocal(uint PosA, uint PosB, const string &Path)
	{
	uint NA = GetNA(Path);
	uint NB = GetNB(Path);
	uint HiA = PosA + NA - 1;
	uint HiB = PosB + NB - 1;
	Log("%3u .. %3u,  %3u .. %3u  %s\n", PosA, HiA, PosB, HiB, Path.c_str());
	}

void cmd_test()
	{
	//EnumPathsGlobal(3, 3, OnPathGlobal);
	EnumPathsLocal(3, 3, OnPathLocal);
	}
