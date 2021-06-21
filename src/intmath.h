// IntMath.h: Header for doing fractional math with integers for speed.

#ifndef IntMath_h
#define	IntMath_h

typedef float BASETYPE;
//typedef double BASETYPE;

// Scaling factor used to store certain floating point
// values as integers to a few significant figures.
//const int INTSCALE = 1000;
const int INTSCALE = 1;

// Type for a probability in range 0.0 to 1.0.
typedef BASETYPE PROB;

// Type for an log-odds integer score.
// Stored as log2(PROB)*INTSCALE.
//typedef int	SCORE;
typedef BASETYPE SCORE;

// Type for a weight.
// Stored as w*INTSCALE where w is in range 0.0 to 1.0.
//typedef unsigned WEIGHT;
typedef BASETYPE WEIGHT;

// Type for a fractional weighted count stored as n*WEIGHT/N
// where n=measured count (integer >= 0) and N is total for
// the distribution (e.g., n=number of residues of a given
// type in a column, N=number of residues in the column).
// Hence values in an FCOUNT variable range from 0..INTSCALE
// as an integer, representing "true" values 0.0 to 1.0.
//typedef unsigned FCOUNT;
typedef BASETYPE FCOUNT;

// Representation of -infinity. Value should
// be large and negative, but not so large
// that adding a few of them overflows.
// TODO: Multiplied by 10 to work around bug
// when aligning Bali 1ckaA in ref4, which is
// so long that B->Mmax got to -infinity, causing
// traceback to fail.
//const int MINUS_INFINITY = -10000000;
const BASETYPE MINUS_INFINITY = (BASETYPE) -1e37;
const BASETYPE PLUS_INFINITY = (BASETYPE) 1e37;

// Probability relative to a null model
typedef double RPROB;

PROB ScoreToProb(SCORE Score);
SCORE ProbToScore(PROB Prob);
SCORE DoubleToScore(double d);
WEIGHT DoubleToWeight(double d);
double WeightToDouble(WEIGHT w);
SCORE MulScoreWeight(SCORE Score, WEIGHT Weight);
bool ScoreEq(SCORE s1, SCORE s2);
bool BTEq(double b1, double b2);

static double ScoreToDouble(SCORE Score)
	{
	return (double) Score / (double) INTSCALE;
	}

#if	0
// In-line assembler for Result = (x*y)/z
// Note that imul and idiv will do 64-bit arithmetic
// on 32-bit operands, so this shouldn't overflow
// Can't write this efficiently in C/C++ (would
// often overlow 32 bits).
#define MulDivAssign(Result, x, y, z)	\
	{									\
	int X = (x);						\
	int Y = (y);						\
	int Z = (z);						\
	_asm mov	eax,X					\
	_asm imul	Y						\
	_asm mov	ecx,Z					\
	_asm idiv	ecx						\
	_asm mov	Result,eax				\
	}
#else
#define MulDivAssign(Result, x, y, z)	Result = (((x)*(y))/(z))
#endif

#define	MulScoreWeight(r, s, w)		MulDivAssign(r, s, w, INTSCALE)
#define MulWeightWCount(r, wt, wc)	MulDivAssign(r, wt, wc, INTSCALE)
#define MulFCountScore(r, fc, sc)	MulDivAssign(r, fc, sc, INTSCALE)

#if	_DEBUG

static inline SCORE Add2(SCORE a, SCORE b)
	{
	if (MINUS_INFINITY == a)
		return MINUS_INFINITY;
	if (MINUS_INFINITY == b)
		return MINUS_INFINITY;
	SCORE sum = a + b;
	if (sum < MINUS_INFINITY)
		return MINUS_INFINITY;
//	assert(sum < OVERFLOW_WARN);
	return sum;
	}

static inline SCORE Add3(SCORE a, SCORE b, SCORE c)
	{
	return Add2(Add2(a, b), c);
	}

static inline SCORE Add4(SCORE a, SCORE b, SCORE c, SCORE d)
	{
	return Add2(Add2(a, b), Add2(c, d));
	}

static inline SCORE Add5(SCORE a, SCORE b, SCORE c, SCORE d, SCORE e)
	{
	return Add3(Add2(a, b), Add2(c, d), e);
	}

static inline SCORE Add6(SCORE a, SCORE b, SCORE c, SCORE d, SCORE e, SCORE f)
	{
	return Add3(Add2(a, b), Add2(c, d), Add2(e, f));
	}

static inline SCORE Add7(SCORE a, SCORE b, SCORE c, SCORE d, SCORE e, SCORE f, SCORE g)
	{
	return Add4(Add2(a, b), Add2(c, d), Add2(e, f), g);
	}

static inline SCORE Mul2(SCORE a, SCORE b)
	{
	if (MINUS_INFINITY == a)
		return MINUS_INFINITY;
	if (MINUS_INFINITY == b)
		return MINUS_INFINITY;
	//__int64 prod = (__int64) a * (__int64) b;
	//assert((SCORE) prod == prod);
	//return (SCORE) prod;
	return a*b;
	}

static inline SCORE Sub2(SCORE a, SCORE b)
	{
	if (MINUS_INFINITY == a)
		return MINUS_INFINITY;
	if (MINUS_INFINITY == b)
		return MINUS_INFINITY;
	SCORE diff = a - b;
	if (diff < MINUS_INFINITY)
		return MINUS_INFINITY;
//	assert(diff < OVERFLOW_WARN);
	return diff;
	}

static inline SCORE Div2(SCORE a, int b)
	{
	if (MINUS_INFINITY == a)
		return MINUS_INFINITY;
	return a/b;
	}

//static inline SCORE MulScoreWeight(SCORE s, WEIGHT w)
//	{
//	SCORE Prod = s*(SCORE) w;
//	assert(Prod < OVERFLOW_WARN);
//	extern void Log(const char Format[], ...);
//	if (Prod/(SCORE) w != s)
//		Log("**WARRNING MulScoreWeight Prod=%d w=%d Prod/w=%d s=%d\n",
//		  Prod,
//		  w,
//		  Prod/(SCORE) w,
//		  s);
//	assert(Prod/ (SCORE) w == s);
//	return Prod/INTSCALE;
//	}
//
//static inline WCOUNT MulWeightWCount(WEIGHT wt, WCOUNT wc)
//	{
//	return (wt*wc)/INTSCALE;
//	}

#else
#define	Add2(a, b)					((a) + (b))
#define Sub2(a, b)					((MINUS_INFINITY == (a)) ? MINUS_INFINITY : ((a) - (b)))
#define Div2(a, b)					((MINUS_INFINITY == (a)) ? MINUS_INFINITY : ((a) / (b)))
#define	Add3(a, b, c)				((a) + (b) + (c))
#define	Add4(a, b, c, d)			((a) + (b) + (c) + (d))
#define	Add5(a, b, c, d, e)			((a) + (b) + (c) + (d) + (e))
#define	Add6(a, b, c, d, e, f)		((a) + (b) + (c) + (d) + (e) + (f))
#define	Add7(a, b, c, d, e, f, g)	((a) + (b) + (c) + (d) + (e) + (f) + (g))
//#define	MulScoreWeight(s, w)		(((s)*(SCORE) (w))/INTSCALE)
#define	Mul2(a, b)					((a)*(b))
#endif

//static inline SCORE MulFCountScore(FCOUNT fc, SCORE sc)
//	{
//// Fast way to say "if (fc >= 2^15 || sc >= 2^15)":
//	if ((fc | sc) & 0xffff1000)
//		{
//		SCORE Score = ((fc+5)/10)*sc;
//		assert(Score < assert);
//		OVERFLOW_WARN(Score > MINUS_INFINITY);
//		return Score/(INTSCALE/10);
//		}
//	SCORE Score = fc*sc;
//	assert(Score < OVERFLOW_WARN);
//	assert(Score > MINUS_INFINITY);
//	return Score/INTSCALE;
//	}

#endif	// IntMath_h
