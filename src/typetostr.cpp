#include "muscle.h"
#include <stdio.h>

const char *SecsToStr(unsigned long Secs)
	{
	static char Str[64];
	long hh, mm, ss;

	hh = Secs/(60*60);
	mm = (Secs/60)%60;
	ss = Secs%60;

	sprintf(Str, "%02ld%c%02ld%c%02ld", hh, ':', mm, ':', ss);
	return Str;
	}

const char *BoolToStr(bool b)
	{
	return b ? "True" : "False";
	}

const char *ScoreToStr(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return "       *";
// Hack to use "circular" buffer so when called multiple
// times in a printf-like argument list it works OK.
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static char szStr[iBufferCount*iBufferLength];
	static int iBufferIndex = 0;
	iBufferIndex = (iBufferIndex + 1)%iBufferCount;
	char *pStr = szStr + iBufferIndex*iBufferLength;
	sprintf(pStr, "%8g", Score);
	return pStr;
	}

// Left-justified version of ScoreToStr
const char *ScoreToStrL(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return "*";
// Hack to use "circular" buffer so when called multiple
// times in a printf-like argument list it works OK.
	const int iBufferCount = 16;
	const int iBufferLength = 16;
	static char szStr[iBufferCount*iBufferLength];
	static int iBufferIndex = 0;
	iBufferIndex = (iBufferIndex + 1)%iBufferCount;
	char *pStr = szStr + iBufferIndex*iBufferLength;
	sprintf(pStr, "%.3g", Score);
	return pStr;
	}

const char *WeightToStr(WEIGHT w)
	{
	return ScoreToStr(w);
	}
