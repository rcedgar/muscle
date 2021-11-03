#include "muscle.h"
#include "ensemble.h"

void cmd_disperse()
	{
	const string FileName = opt(disperse);

	Ensemble E;
	E.FromFile(FileName);

	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	double D_LP;
	double D_Cols;
	E.GetDispersion(MaxGapFract, D_LP, D_Cols);

	ProgressLog("@disperse file=%s D_LP=%.4g D_Cols=%.4g\n",
	  FileName.c_str(), D_LP, D_Cols);
	}
