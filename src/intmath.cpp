#include "muscle.h"
#include <math.h>

PROB ScoreToProb(SCORE Score)
	{
	if (MINUS_INFINITY >= Score)
		return 0.0;
	return (PROB) pow(2.0, (double) Score/INTSCALE);
	}

//#if	0
//static const double log2e = log2(exp(1.0));
//
//double lnTolog2(double ln)
//	{
//	return ln*log2e;
//	}
//
//double log2(double x)
//	{
//	if (0 == x)
//		return MINUS_INFINITY;
//
//	static const double dInvLn2 = 1.0/log(2.0);
//// Multiply by inverse of log(2) just in case multiplication
//// is faster than division.
//	return log(x)*dInvLn2;
//	}
//#endif

//SCORE ProbToScore(PROB Prob)
//	{
//	if (0.0 == Prob)
//		return MINUS_INFINITY;
////	return (SCORE) floor(INTSCALE*log2(Prob));
//	return (SCORE) log2(Prob);
//	}

WEIGHT DoubleToWeight(double d)
	{
	assert(d >= 0);
	return (WEIGHT) (INTSCALE*d);
	}

double WeightToDouble(WEIGHT w)
	{
	return (double) w / (double) INTSCALE;
	}

SCORE DoubleToScore(double d)
	{
	return (SCORE)(d*(double) INTSCALE);
	}

bool ScoreEq(SCORE s1, SCORE s2)
	{
	return BTEq(s1, s2);
	}

static bool BTEq2(BASETYPE b1, BASETYPE b2)
	{
	double diff = fabs(b1 - b2);
	if (diff < 0.0001)
		return true;
	double sum = fabs(b1) + fabs(b2);
	return diff/sum < 0.005;
	}

bool BTEq(double b1, double b2)
	{
	return BTEq2((BASETYPE) b1, (BASETYPE) b2);
	}

//const double dLn2 = log(2.0);

//// pow2(x)=2^x
//double pow2(double x)
//	{
//	if (MINUS_INFINITY == x)
//		return 0;
//	return exp(x*dLn2);
//	}

//// lp2(x) = log2(1 + 2^-x), x >= 0
//double lp2(double x)
//	{
//	return log2(1 + pow2(-x));
//	}

// SumLog(x, y) = log2(2^x + 2^y)
//SCORE SumLog(SCORE x, SCORE y)
//	{
//	return (SCORE) log2(pow2(x) + pow2(y));
//	}
//
//// SumLog(x, y, z) = log2(2^x + 2^y + 2^z)
//SCORE SumLog(SCORE x, SCORE y, SCORE z)
//	{
//	return (SCORE) log2(pow2(x) + pow2(y) + pow2(z));
//	}
//
//// SumLog(w, x, y, z) = log2(2^w + 2^x + 2^y + 2^z)
//SCORE SumLog(SCORE w, SCORE x, SCORE y, SCORE z)
//	{
//	return (SCORE) log2(pow2(w) + pow2(x) + pow2(y) + pow2(z));
//	}

//SCORE lp2Fast(SCORE x)
//	{
//	assert(x >= 0);
//	const int iTableSize = 1000;
//	const double dRange = 20.0;
//	const double dScale = dRange/iTableSize;
//	static SCORE dValue[iTableSize];
//	static bool bInit = false;
//	if (!bInit)
//		{
//		for (int i = 0; i < iTableSize; ++i)
//			dValue[i] = (SCORE) lp2(i*dScale);
//		bInit = true;
//		}
//	if (x >= dRange)
//		return 0.0;
//	int i = (int) (x/dScale);
//	assert(i >= 0 && i < iTableSize);
//	SCORE dResult = dValue[i];
//	assert(BTEq(dResult, lp2(x)));
//	return dResult;
//	}
//
//// SumLog(x, y) = log2(2^x + 2^y)
//SCORE SumLogFast(SCORE x, SCORE y)
//	{
//	if (MINUS_INFINITY == x)
//		{
//		if (MINUS_INFINITY == y)
//			return MINUS_INFINITY;
//		return y;
//		}
//	else if (MINUS_INFINITY == y)
//		return x;
//
//	SCORE dResult;
//	if (x > y)
//		dResult = x + lp2Fast(x-y);
//	else
//		dResult = y + lp2Fast(y-x);
//	assert(SumLog(x, y) == dResult);
//	return dResult;
//	}
//
//SCORE SumLogFast(SCORE x, SCORE y, SCORE z)
//	{
//	SCORE dResult = SumLogFast(x, SumLogFast(y, z));
//	assert(SumLog(x, y, z) == dResult);
//	return dResult;
//	}

//SCORE SumLogFast(SCORE w, SCORE x, SCORE y, SCORE z)
//	{
//	SCORE dResult = SumLogFast(SumLogFast(w, x), SumLogFast(y, z));
//	assert(SumLog(w, x, y, z) == dResult);
//	return dResult;
//	}

double VecSum(const double v[], unsigned n)
	{
	double dSum = 0.0;
	for (unsigned i = 0; i < n; ++i)
		dSum += v[i];
	return dSum;
	}

void Normalize(PROB p[], unsigned n)
	{
	unsigned i;
	PROB dSum = 0.0;
	for (i = 0; i < n; ++i)
		dSum += p[i];
	if (0.0 == dSum)
		Quit("Normalize, sum=0");
	for (i = 0; i < n; ++i)
		p[i] /= dSum;
	}

void NormalizeUnlessZero(PROB p[], unsigned n)
	{
	unsigned i;
	PROB dSum = 0.0;
	for (i = 0; i < n; ++i)
		dSum += p[i];
	if (0.0 == dSum)
		return;
	for (i = 0; i < n; ++i)
		p[i] /= dSum;
	}

void Normalize(PROB p[], unsigned n, double dRequiredTotal)
	{
	unsigned i;
	double dSum = 0.0;
	for (i = 0; i < n; ++i)
		dSum += p[i];
	if (0.0 == dSum)
		Quit("Normalize, sum=0");
	double dFactor = dRequiredTotal / dSum;
	for (i = 0; i < n; ++i)
		p[i] *= (PROB) dFactor;
	}

bool VectorIsZero(const double dValues[], unsigned n)
	{
	for (unsigned i = 0; i < n; ++i)
		if (dValues[i] != 0.0)
			return false;
	return true;
	}

void VectorSet(double dValues[], unsigned n, double d)
	{
	for (unsigned i = 0; i < n; ++i)
		dValues[i] = d;
	}

bool VectorIsZero(const float dValues[], unsigned n)
	{
	for (unsigned i = 0; i < n; ++i)
		if (dValues[i] != 0.0)
			return false;
	return true;
	}

void VectorSet(float dValues[], unsigned n, float d)
	{
	for (unsigned i = 0; i < n; ++i)
		dValues[i] = d;
	}

double Correl(const double P[], const double Q[], unsigned uCount)
	{
	double dSumP = 0.0;
	double dSumQ = 0.0;
	for (unsigned n = 0; n < uCount; ++n)
		{
		dSumP += P[n];
		dSumQ += Q[n];
		}
	const double dMeanP = dSumP/uCount;
	const double dMeanQ = dSumQ/uCount;

	double dSum1 = 0.0;
	double dSum2 = 0.0;
	double dSum3 = 0.0;
	for (unsigned n = 0; n < uCount; ++n)
		{
		const double dDiffP = P[n] - dMeanP;
		const double dDiffQ = Q[n] - dMeanQ;
		dSum1 += dDiffP*dDiffQ;
		dSum2 += dDiffP*dDiffP;
		dSum3 += dDiffQ*dDiffQ;
		}
	if (0 == dSum1)
		return 0;
	const double dCorrel = dSum1 / sqrt(dSum2*dSum3);
	return dCorrel;
	}

float Correl(const float P[], const float Q[], unsigned uCount)
	{
	float dSumP = 0.0;
	float dSumQ = 0.0;
	for (unsigned n = 0; n < uCount; ++n)
		{
		dSumP += P[n];
		dSumQ += Q[n];
		}
	const float dMeanP = dSumP/uCount;
	const float dMeanQ = dSumQ/uCount;

	float dSum1 = 0.0;
	float dSum2 = 0.0;
	float dSum3 = 0.0;
	for (unsigned n = 0; n < uCount; ++n)
		{
		const float dDiffP = P[n] - dMeanP;
		const float dDiffQ = Q[n] - dMeanQ;
		dSum1 += dDiffP*dDiffQ;
		dSum2 += dDiffP*dDiffP;
		dSum3 += dDiffQ*dDiffQ;
		}
	if (0 == dSum1)
		return 0;
	const float dCorrel = dSum1 / (float) sqrt(dSum2*dSum3);
	return dCorrel;
	}

// Simple (but slow) function to compute Pearson ranks
// that allows for ties. Correctness and simplicity
// are priorities over speed here.
void Rank(const float P[], float Ranks[], unsigned uCount)
	{
	for (unsigned n = 0; n < uCount; ++n)
		{
		unsigned uNumberGreater = 0;
		unsigned uNumberEqual = 0;
		unsigned uNumberLess = 0;
		double dValue = P[n];
		for (unsigned i = 0; i < uCount; ++i)
			{
			double v = P[i];
			if (v == dValue)
				++uNumberEqual;
			else if (v < dValue)
				++uNumberLess;
			else
				++uNumberGreater;
			}
		assert(uNumberEqual >= 1);
		assert(uNumberEqual + uNumberLess + uNumberGreater == uCount);
		Ranks[n] = (float) (1 + uNumberLess + (uNumberEqual - 1)/2.0);
		}
	}

void Rank(const double P[], double Ranks[], unsigned uCount)
	{
	for (unsigned n = 0; n < uCount; ++n)
		{
		unsigned uNumberGreater = 0;
		unsigned uNumberEqual = 0;
		unsigned uNumberLess = 0;
		double dValue = P[n];
		for (unsigned i = 0; i < uCount; ++i)
			{
			double v = P[i];
			if (v == dValue)
				++uNumberEqual;
			else if (v < dValue)
				++uNumberLess;
			else
				++uNumberGreater;
			}
		assert(uNumberEqual >= 1);
		assert(uNumberEqual + uNumberLess + uNumberGreater == uCount);
		Ranks[n] = (double) (1 + uNumberLess + (uNumberEqual - 1)/2.0);
		}
	}

FCOUNT SumCounts(const FCOUNT Counts[])
	{
	FCOUNT Sum = 0;
	for (int i = 0; i < 20; ++i)
		Sum += Counts[i];
	return Sum;
	}
