#include "muscle.h"

uint64 GetFBSize(uint LX, uint LY)
	{
	uint64 Size64 = uint64(LX + 1)*uint64(LY + 1)*HMMSTATE_COUNT;
	return Size64;
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
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		Die("Sequences length %u, %u overflow HMM buffers", LX, LY);

	uint64 Bytes = GetFBSize(LX, LY)*uint64(sizeof(float));
	return myalloc(float, Bytes);
	}

float *AllocPost(uint LX, uint LY)
	{
	return myalloc(float, GetPostSize(LX, LY));
	}

float *AllocDPRows(uint LX, uint LY)
	{
	return myalloc(float, GetDPRowsSize(LX, LY));
	}

char *AllocTB(uint LX, uint LY)
	{
	return myalloc(char, GetTBSize(LX, LY));
	}
