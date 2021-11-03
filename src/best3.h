#pragma once

// Store the largest of three values x1, x2, and x3 in *x.
// If x_i is the largest value, then store b_i in *b.
static inline void Best3(float x1, float x2, float x3,
  char b1, char b2, char b3, float *x, char *b)
	{
	if (x1 >= x2)
		{
		if (x1 >= x3)
			{
			*x = x1;
			*b = b1;
			return;
			}
		*x = x3;
		*b = b3;
		return;
		}
	if (x2 >= x3)
		{
		*x = x2;
		*b = b2;
		return;
		}
	*x = x3;
	*b = b3;
	}

// Store the largest of three values x1, x2, and x3 in *x.
static inline void Best3(float x1, float x2, float x3, float *x)
	{
	if (x1 >= x2)
		{
		if (x1 >= x3)
			{
			*x = x1;
			return;
			}
		*x = x3;
		return;
		}
	if (x2 >= x3)
		{
		*x = x2;
		return;
		}
	*x = x3;
	}
