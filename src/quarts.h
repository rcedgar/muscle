#ifndef quarts_h
#define quarts_h

#include "myutils.h"

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
	unsigned N;
	float Min;
	float LoQ;
	float Med;
	float HiQ;
	float Max;
	float Total;
	float Avg;
	float StdDev;

	void ProgressLogMe() const
		{
		ProgressLog("N=%u", N);
		ProgressLog(", Min=%.3g", Min);
		ProgressLog(", LoQ=%.3g", LoQ);
		ProgressLog(", Med=%.3g", Med);
		ProgressLog(", HiQ=%.3g", HiQ);
		ProgressLog(", Max=%.3g", Max);
		ProgressLog(", Avg=%.3g", Avg);
		ProgressLog(", StdDev=%.3g", StdDev);
		ProgressLog("\n");
		}

	void ToTsv(FILE *f) const
		{
		if (f == 0)
			return;
		fprintf(f, "%u", N);
		fprintf(f, "\t%.3g", Min);
		fprintf(f, "\t%.3g", LoQ);
		fprintf(f, "\t%.3g", Med);
		fprintf(f, "\t%.3g", HiQ);
		fprintf(f, "\t%.3g", Max);
		fprintf(f, "\t%.3g", Avg);
		fprintf(f, "\t%.3g", StdDev);
		fprintf(f, "\n");
		}
	};

void GetQuarts(const vector<unsigned> &v, Quarts &Q);
void GetQuartsFloat(const vector<float> &v, QuartsFloat &Q);

#endif // quarts_h
