// scoretype.h
#pragma once

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

/***
https://gasstationwithoutpumps.wordpress.com/2014/05/06/sum-of-probabilities-in-log-prob-space/

The problem is mathematically very simple.  We have two non-negative numbers
(usually probabilities), p and q, that are represented in the computer by
floating-point numbers for their natural logs: a=\ln p and b=\ln q, and we want
to represent their sum: \ln(p+q).

The brute-force way to do this in code is log(exp(a) + exp(b)), but this can run
into floating-point problems (underflow, overflow, or loss of precision). It is
also fairly expensive, as it involves computing three transcendental functions
(2 exponentiation and 1 logarithm).

We can simplify the problem:
	\ln(e^{a} + e^{b})
		= \ln\left( e^{b}(e^{a-b} + 1) \right)
		= b + \ln(e^{a-b} +1),

so all we need to compute is the function f(x) = \ln(e^{x} +1).  Furthermore,
if we swap the inputs as needed, we can make sure that b\geq a, so that f(x)
is only needed for x\leq 0.  This means that we only need two transcendental
functions (one exponential and one logarithm), and the logarithm is always
going to be of a number between 1 and 2 (which is generally where logarithm
implementations are most efficient). We can get much better accuracy by
eliminating the addition and using library function log1p which computes
\ln(1+x). Our function f is then computed as log1p(exp(x)).

But we can still run into problems.  If x is very negative, then exp(x)
could underflow to 0. This is not a serious problem, unless the programming
language treats that as an exception and throws an error. We can avoid this
problem by setting a threshold, below which the implementation of f just
returns 0. The threshold should be set at the most negative value of x
for which exp(x) still returns a normal floating-point number: that is,
the natural log of the smallest normalized number, which depends on the
precision. For the IEEE standard floating-point representations, I believe
that the smallest normalized numbers are

representation 	smallest +vs float 	approx ln
float16 	    2^-15 	            -10.397207708
float32 	    2^-127 	            -88.029691931
float64 	    2^-1023 	        -709.089565713
float128 	    2^-16383 	        -11355.8302591

By adding one test for the threshold (appropriately chosen for the precision
of floating-point number used in the log-prob representation), and returning
0 from f when below the threshold, we can avoid underflowing on the computation
of the exponential.

I was going to write up the choice of a slightly higher cutpoint, where we could
use the Taylor expansion \ln(1+\epsilon) \approx \epsilon + O(\epsilon^2), and
just return exp(x), but I believe that log1p already handles this correctly.
So reasonably efficient code that does the sum of probabilities in log-prob 
representation looks something like this in c++:

inline float log1pexp(float x)
{   return x<-88.029691931? 0.: log1p(exp(x));
}
inline float sum_log_prob(float a, float b)
{   return a>b? a+log1pexp(b-a):  b+log1pexp(a-b);
}
inline double log1pexp(double x)
{   return x<-709.089565713? 0.: log1p(exp(x));
}
inline double sum_log_prob(double a, double b)
{   return a>b? a+log1pexp(b-a):  b+log1pexp(a-b);
}
inline long double log1pexp(long double x)
{   return x<-11355.8302591? 0.: log1p(exp(x));
}
inline long double sum_log_prob(long double a, long double b)
{   return a>b? a+log1pexp(b-a):  b+log1pexp(a-b);
}

A careful coder would check exactly where the exp(x) computation underflows,
rather than relying on the theoretical estimates made here, as there could
be implementation-dependent details that cause underflow at a higher threshold
than strictly necessary.
***/

const float LOG_ZERO = -2e20f;
const float LOG_ONE = 0.0f;
const float INVALID_LOG = FLT_MAX;
const float UNINIT_LOG = 9e9f;
const float OUT_OF_BAND_LOG = 8e8f;

const float EXP_UNDERFLOW_THRESHOLD = -4.6f;
const float LOG_UNDERFLOW_THRESHOLD = 7.5f;

// Computes log (exp (x) + 1), for 0 <= x <= 7.5.
// This is 2x faster than the libary function log1p()
inline float LOGEXP1(float x)
	{
	assert(x >= 0.00f);
	assert(x <= LOG_UNDERFLOW_THRESHOLD);
	if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
	if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
	if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
	assert(x <= LOG_UNDERFLOW_THRESHOLD);
	return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
	}

inline void LOG_PLUS_EQUALS(float& x, float y)
	{
	if (x < y)
		x = (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOGEXP1(y - x) + x;
	else
		x = (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOGEXP1(x - y) + y;
	}

inline float LOG_ADD(float x, float y)
	{
	if (x < y)
		return (x == LOG_ZERO || y - x >= LOG_UNDERFLOW_THRESHOLD) ? y : LOGEXP1(y - x) + x;
	return (y == LOG_ZERO || x - y >= LOG_UNDERFLOW_THRESHOLD) ? x : LOGEXP1(x - y) + y;
	}

inline float LOG_ADD(float x1, float x2, float x3)
	{
	return LOG_ADD(x1, LOG_ADD(x2, x3));
	}

inline float LOG_ADD(float x1, float x2, float x3, float x4)
	{
	return LOG_ADD(x1, LOG_ADD(x2, LOG_ADD(x3, x4)));
	}

inline float LOG_ADD(float x1, float x2, float x3, float x4, float x5)
	{
	return LOG_ADD(x1, LOG_ADD(x2, LOG_ADD(x3, LOG_ADD(x4, x5))));
	}

inline float LOG_ADD(float x1, float x2, float x3, float x4, float x5, float x6)
	{
	return LOG_ADD(x1, LOG_ADD(x2, LOG_ADD(x3, LOG_ADD(x4, LOG_ADD(x5, x6)))));
	}

inline float LOG_ADD(float x1, float x2, float x3, float x4, float x5, float x6, float x7)
	{
	return LOG_ADD(x1, LOG_ADD(x2, LOG_ADD(x3, LOG_ADD(x4, LOG_ADD(x5, LOG_ADD(x6, x7))))));
	}
