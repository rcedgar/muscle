#include "myutils.h"
#include "alpha.h"

/***
This source file adapted from PHYLIP 3.696:
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Dan Fineman.

   Copyright (c) 1993-2014, Joseph Felsenstein
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
   POSSIBILITY OF SUCH DAMAGE.
***/

void RePredict(uint Letterq, uint Lettert,
  double *tt, double *p, double *dp, double *d2p, double *q, double *elambdat);

const uint MIN_OVERLAP = 8;
const uint ITERS = 20;
const double protepsilon = 0.00001;

double GetProtDist(const char *Q, const char *T, uint ColCount)
	{
	double lnlike, slope, curv;
	bool neginfinity, inf, overlap;
	double p, dp, d2p;
	double elambdat = 0;
	double q = 0;

	double tt = 0.1;
	double delta = tt / 2.0;
	inf = false;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		lnlike = 0.0;
		slope = 0.0;
		curv = 0.0;
		neginfinity = false;
		overlap = false;
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			byte CharQ = (byte) Q[Col];
			byte CharT = (byte) T[Col];
			uint LetterQ = g_CharToLetterAmino[CharQ];
			uint LetterT = g_CharToLetterAmino[CharT];
			if (LetterQ >= 20 || LetterT >= 20)
				continue;
			overlap = true;
			p = 0.0;
			dp = 0.0;
			d2p = 0.0;
			RePredict(LetterQ, LetterT, &tt, &p, &dp, &d2p, &q, &elambdat);
			if (p <= 0.0)
				neginfinity = true;
			else
				{
				lnlike += log(p);
				slope += dp / p;
				curv += (d2p / p - dp * dp / (p * p));
				}
			} // end for Col
		if (!overlap)
			return -1;
		else if (!neginfinity)
			{
			if (curv < 0.0)
				{
				tt -= slope / curv;
				if (tt > 10000.0)
					return -1;
				}
			else
				{
				if ((slope > 0.0 && delta < 0.0) || (slope < 0.0 && delta > 0.0))
				delta /= -2;
				tt += delta;
				}
			}
		else
			{
			delta /= -2;
			tt += delta;
			}
		if (tt < protepsilon && !inf)
			tt = protepsilon;
		} // end for Iter
	double d = tt;
	return d;
	}

double GetProtDist(const string &Q, const string &T)
	{
	uint ColCount = SIZE(Q);
	asserta(SIZE(T) == ColCount);
	double d = GetProtDist(Q.c_str(), T.c_str(), ColCount);
	return d;
	}
