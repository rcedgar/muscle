#include "muscle.h"

// BLASTZ default parameters
// open 400, extend 30, matrix as below

const float NUC_EXTEND = 30;
const float NUC_SP_CENTER = 2*NUC_EXTEND;

#define v(x)	((float) x + NUC_SP_CENTER)
#define ROW(A, C, G, T) \
	{ v(A), v(C), v(G), v(T) },

float NUC_SP[32][32] =
	{
//         A        C        G        T
ROW(      91,    -114,     -31,    -123) // A

ROW(    -114,     100,    -125,     -31) // C

ROW(     -31,    -125,     100,    -114) // G

ROW(    -123,     -31,    -114,      91) // T
	};
