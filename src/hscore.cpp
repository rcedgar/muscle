#include "myutils.h"

double hscore(const double *xs, const double *ys, uint N, double X)
	{
	if (X < xs[0])
		return ys[0];

	for (uint i = 1; i < N; ++i)
		{
		if (X >= xs[i-1] && X < xs[i])
			{
			double f = (X - xs[i-1])/(xs[i] - xs[i-1]);
			double y = f*ys[i] + (1 - f)*ys[i-1];
			return y;
			}
		}
	return ys[N-1];
	}

#if 0
static void Test(const double *xs, const double *ys, uint N, double X)
	{
	double y = hscore(xs, ys, N, X);
	Log("x=%.3g y=%.3g\n", X, y);
	}

void cmd_test()
	{
	const double xs[] = {0.1, 0.5, 0.9};
	const double ys[] = {1, 2, 3};
	const uint N = sizeof(xs)/sizeof(xs[0]);
	asserta(sizeof(xs)/sizeof(xs[0]) == N);

	Test(xs, ys, N, 0.0);
	Test(xs, ys, N, 0.1);
	Test(xs, ys, N, 0.2);
	Test(xs, ys, N, 0.3);
	Test(xs, ys, N, 0.4);
	Test(xs, ys, N, 0.5);
	Test(xs, ys, N, 0.6);
	Test(xs, ys, N, 0.7);
	Test(xs, ys, N, 0.8);
	Test(xs, ys, N, 0.9);
	Test(xs, ys, N, 1.0);
	}
#endif // 0
