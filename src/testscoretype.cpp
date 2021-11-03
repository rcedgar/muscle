#include "muscle.h"
#include "timing.h"

void Test2(double P1d, double P2d)
	{
	float P1 = float(P1d);
	float P2 = float(P2d);

	float Prod = P1*P2;
	float Sum = P1 + P2;

	float logP1 = log(P1);
	float logP2 = log(P2);
	float logProd = log(Prod);
	float logSum = log(Sum);

	float sumP12 = logP1 + logP2;
	asserta(feq(logProd, sumP12));

	float Add = LOG_ADD(logP1, logP2);
	asserta(feq(Add, logSum));

	float PE = logP1;
	LOG_PLUS_EQUALS(PE, logP2);
	asserta(feq(PE, logSum));

	ProgressLog("P1=%.3g P2=%.3g ok\n", P1, P2);
	}

#if 0
static void TestExp()
	{
	vector<float> Xs;
	Progress("Calc Xs...");
	for (uint i = 0; i < 1000000; ++i)
		{
		uint r = randu32()%16;
		float x = -float(r)/16;
		Xs.push_back(x);
		}
	Progress("\n");

	TICKS t1 = GetClockTicks();
	float Sum = 0;
	for (uint i = 0; i < SIZE(Xs); ++i)
		{
		float x = Xs[i];
		float E = EXP(x);
		Sum += E;
		}
	TICKS t2 = GetClockTicks();

	TICKS t3 = GetClockTicks();
	for (uint i = 0; i < SIZE(Xs); ++i)
		{
		float x = Xs[i];
		float e = exp(x);
		Sum += e;
		}
	TICKS t4 = GetClockTicks();

	ProgressLog("EXP %.3g, exp %.3g\n", double(t2 - t1), double(t4 - t3));
	}
#endif // 0
