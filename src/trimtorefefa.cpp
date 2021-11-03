#include "myutils.h"
#include "alpha.h"
#include "msa.h"
#include "ensemble.h"

void TrimToRef(const MSA &Test, const MSA &Ref, MSA &Trimmed);

void cmd_trimtoref_efa()
	{
	const string EfaFileName = opt(trimtoref_efa);
	const string RefFileName = opt(ref);
	const string OutputFileName = opt(output);

	Ensemble E;
	E.FromFile(EfaFileName);

	MSA Ref;
	Ref.FromFASTAFile_PreserveCase(RefFileName);

	FILE *fOut = CreateStdioFile(OutputFileName);

	const uint MSACount = E.GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		Pf(fOut, "<%s\n", E.GetMSAName(MSAIndex).c_str());
		const MSA &TestMSA = E.GetMSA(MSAIndex);
		MSA TrimmedMSA;
		TrimToRef(TestMSA, Ref, TrimmedMSA);
		TrimmedMSA.ToFASTAFile(fOut);
		}

	CloseStdioFile(fOut);
	}
