#include "muscle.h"

void LogTomMx(const string &Name, const vector<float> &Mx, uint LX, uint LY)
	{
	Log("\n");
	Log("Tom %s: LX=%u LY=%u\n", Name.c_str(), LX, LY);
	Log("       ");
	for (uint j = 0; j <= LY; ++j)
		Log("  %10u", j);
	Log("\n");

	uint Ix = 0;
	for (uint i = 0; i <= LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j <= LY; ++j)
			{
			float P = Mx[Ix++];
			Log("  %10.3g", P);
			}
		Log("\n");
		}
	}

// (LX+1) x (LY+1)
void LogFlatMx1(const string &Name, const float *MyPost, uint LX, uint LY)
	{
	Log("\n");
	Log("Flat1 %s: LX=%u LY=%u\n", Name.c_str(), LX, LY);
	Log("       ");
	for (uint j = 0; j <= LY; ++j)
		Log("  %10u", j);
	Log("\n");

	uint Ix = 0;
	for (uint i = 0; i <= LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j <= LY; ++j)
			{
			float P = MyPost[Ix++];
			Log("  %10.3g", P);
			}
		Log("\n");
		}
	}

// LX x LY
void LogFlatMx(const string &Name, const float *MyPost, uint LX, uint LY)
	{
	Log("\n");
	Log("Flat %s: LX=%u LY=%u\n", Name.c_str(), LX, LY);
	Log("       ");
	for (uint j = 0; j < LY; ++j)
		Log("  %10u", j);
	Log("\n");

	uint Ix = 0;
	for (uint i = 0; i < LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j < LY; ++j)
			{
			float P = MyPost[Ix++];
			Log("  %10.3g", P);
			}
		Log("\n");
		}
	}

// 5 x (LX + 1) x (LY + 1)
void LogFlatMxs(const string &Name, const float *Mxs, uint LX, uint LY)
	{
	Log("\n");

	for (uint s = 0; s < HMMSTATE_COUNT; ++s)
		{
		Log("Flat %s[%u]: LX=%u LY=%u\n", Name.c_str(), s, LX, LY);
		Log("       ");
		for (uint j = 0; j <= LY; ++j)
			Log("  %10u", j);
		Log("\n");

		uint Ix = s;
		for (uint i = 0; i <= LX; ++i)
			{
			Log("[%3u]  ", i);
			for (uint j = 0; j <= LY; ++j)
				{
				float x = Mxs[Ix];
				Ix += HMMSTATE_COUNT;
				if (x == INVALID_LOG)
					Log("  %8.8s", "*ERR*");
				if (x == OUT_OF_BAND_LOG)
					Log("  %8.8s", "#");
				if (x == UNINIT_LOG)
					Log("  %8.8s", "-");
				else if (x == LOG_ZERO)
					Log("  %8.8s", ".");
				else
					Log("  %8.3g", x);
				}
			Log("\n");
			}
		}
	}
