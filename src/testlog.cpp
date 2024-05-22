#include "muscle.h"
#include "timing.h"

//void cmd_testlog() {}

#if 0

inline float HACK(float x)
	{
	assert(x >= 0.00f);
	assert(x <= LOG_UNDERFLOW_THRESHOLD);
	if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
	//if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
	//if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
	assert(x <= LOG_UNDERFLOW_THRESHOLD);
	return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
	}

inline float LOG_ADD_HACK(float x, float y)
	{
	if (x < y)
		return (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : HACK(y - x) + x;
	return (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : HACK(x - y) + y;
	}

inline float log1pexp(float x)
	{
	return x < -88.029691931f? 0.0f : log1p(exp(x));
	}

inline float sum_log_prob(float x, float y)
	{
	return x > y ? x + log1pexp(y - x) :  y + log1pexp(x - y);
	}

void cmd_testlog()
	{
	opt(testlog);

	const uint N = 10000;
	float *v = myalloc(float, N);
	for (uint i = 0; i < N; ++i)
		{
		float x = -float(randu32()%20) + float(randu32()%100)/10000.0f;
		v[i] = x;
		}

	const uint M = 1000*1000;
	float *mv = myalloc(float, M);
	ProgressLog("Making mv...");
	for (uint i = 0; i < M; ++i)
		mv[i] = -float(randu32()%20) + float(randu32()%100)/10000.0f;
	ProgressLog("done.\n");

	uint Diffs = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum1 = LOG_ADD(x, y);
			float Sum2 = sum_log_prob(x, y);
			if (!feq(Sum1, Sum2))
				{
				ProgressLog("x=%.5g y=%.5g LOG_ADD=%.3g sum_log_prob=%.3g\n",
				  x, y, Sum1, Sum2);
				if (++Diffs > 20)
					goto Next;
				}
			}
		}

Next:
	TICKS t1 = GetClockTicks();
	float Total = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum = LOG_ADD(x, y);
			Total += Sum;
			}
		}
	TICKS t2 = GetClockTicks();
	ProgressLog("LOG_ADD %.4g ticks\n", double(t2 - t1));
	Log("%.3g\n", Total);

	TICKS t3 = GetClockTicks();
	Total = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum = sum_log_prob(x, y);
			Total += Sum;
			}
		}
	TICKS t4 = GetClockTicks();
	ProgressLog("sum_log_prob %.4g ticks\n", double(t4 - t3));
	Log("%.3g\n", Total);

	TICKS t5 = GetClockTicks();
	Total = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum = mv[randu32()%M];
			Total += Sum;
			}
		}
	Log("%.3g\n", Total);
	TICKS t6 = GetClockTicks();
	ProgressLog("lookup %.4g ticks\n", double(t6 - t5));

	TICKS t7 = GetClockTicks();
	Total = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum = LOG_ADD_HACK(x, y);
			Total += Sum;
			}
		}
	Log("%.3g\n", Total);
	TICKS t8 = GetClockTicks();
	ProgressLog("HACK %.4g ticks\n", double(t8 - t7));

	TICKS t9 = GetClockTicks();
	Total = 0;
	for (uint i = 0; i < N; ++i)
		{
		for (uint j = 0; j < N; ++j)
			{
			float x = v[i];
			float y = v[j];
			float Sum = x + y;
			Total += Sum;
			}
		}
	Log("%.3g\n", Total);
	TICKS t10 = GetClockTicks();
	ProgressLog("NULL %.4g ticks\n", double(t10 - t9));

	ProgressLog("LOG_ADD %.3g sum_log_prob %.3g lookup %.3g\n",
	  double(t2 - t1), double(t4 - t3), double(t6 - t5));
	}
#endif // 0
