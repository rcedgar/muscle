#include "muscle.h"
#include "ensemble.h"

void cmd_efa_explode()
	{
	const string &InputFileName = opt(efa_explode);

	string Prefix;
	if (optset_prefix)
		Prefix = opt(prefix);

	string Suffix;
	if (optset_suffix)
		Suffix = opt(suffix);

	Ensemble E;
	E.FromFile(InputFileName);

	const uint MSACount = E.GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const MSA &M = E.GetMSA(MSAIndex);

		string FileName = E.GetMSAName(MSAIndex);
		if (FileName == "")
			Ps(FileName, "%u", MSAIndex);
		FileName = Prefix + FileName + Suffix;
		M.ToFASTAFile(FileName);
		}
	}
