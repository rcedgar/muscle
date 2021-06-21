#ifndef quarts_h
#define quarts_h

struct Quarts
	{
	unsigned Min;
	unsigned LoQ;
	unsigned Med;
	unsigned HiQ;
	unsigned Max;
	unsigned Total;
	double Avg;
	};

struct QuartsFloat
	{
	float Min;
	float LoQ;
	float Med;
	float HiQ;
	float Max;
	float Total;
	float Avg;
	float StdDev;
	};

void GetQuarts(const vector<unsigned> &v, Quarts &Q);
void GetQuartsFloat(const vector<float> &v, QuartsFloat &Q);

#endif // quarts_h
