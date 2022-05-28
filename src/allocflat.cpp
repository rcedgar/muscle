#include "muscle.h"

uint64 GetFBSize(uint LX, uint LY)
	{
	uint64 Size64 = uint64(LX + 1)*uint64(LY + 1)*HMMSTATE_COUNT;
	if (double(Size64) > 4e9)
		Die("Memory object too large due to sequence lengths %u, %u", LX, LY);
	uint Size = uint(Size64);
	asserta(Size == uint(Size64));
	return Size;
	}

uint64 GetPostSize(uint LX, uint LY)
	{
	uint64 Size64 = uint64(LX)*uint64(LY);
	uint Size = uint(Size64);
	asserta(uint64(Size) == Size64);
	return Size;
	}

uint64 GetDPRowsSize(uint LX, uint LY)
	{
	uint64 Size64 = 2*uint64(LY + 1);
	uint Size = uint(Size64);
	asserta(uint64(Size) == Size64);
	return Size;
	}

uint64 GetTBSize(uint LX, uint LY)
	{
	uint64 Size64 = uint64(LX + 1)*uint64(LY + 1);
	uint Size = uint(Size64);
	asserta(uint64(Size) == Size64);
	return Size;
	}

float *AllocFB(uint LX, uint LY)
	{
	return myalloc64(float, GetFBSize(LX, LY));
	}

float *AllocPost(uint LX, uint LY)
	{
	return myalloc64(float, GetPostSize(LX, LY));
	}

float *AllocDPRows(uint LX, uint LY)
	{
	return myalloc64(float, GetDPRowsSize(LX, LY));
	}

char *AllocTB(uint LX, uint LY)
	{
	return myalloc64(char, GetTBSize(LX, LY));
	}
