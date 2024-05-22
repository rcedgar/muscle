#include "muscle.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "omplock.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI);
void MakeAlnRows(const byte *XSeq, uint LX,
  const byte *YSeq, uint LY, const PathInfo &PI,
  string &RowX, string &RowY);
double GetProtDist(const char *Q, const char *T, uint ColCount);
XDPMem &GetDPMem();

void cmd_protdists()
	{
	const string &InputFileName = opt(protdists);
	const string &OutputFileName = opt(output);
	FILE *fOut = CreateStdioFile(OutputFileName);

	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName);
	bool IsNucleo = InputSeqs.GuessIsNucleo();
	SetAlphab(IsNucleo);

	uint SeqCount = InputSeqs.GetSeqCount();
	uint ThreadCount = GetRequestedThreadCount();
	const uint PairCount = (SeqCount*(SeqCount-1))/2;
	uint SeqIndexi = UINT_MAX;
	uint SeqIndexj = UINT_MAX;
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		uint ThreadLocali = UINT_MAX;
		uint ThreadLocalj = UINT_MAX;

		LOCK();
		if (SeqIndexi == UINT_MAX)
			{
			SeqIndexi = 1;
			SeqIndexj = 0;
			}
		else
			{
			++SeqIndexj;
			if (SeqIndexj == SeqIndexi)
				{
				++SeqIndexi;
				SeqIndexj = 0;
				}
			}
		ThreadLocali = SeqIndexi;
		ThreadLocalj = SeqIndexj;
		ProgressStep(PairCounter++, PairCount, "Protdists");
		UNLOCK();

		uint Li;
		const byte *Seqi = InputSeqs.GetByteSeq(ThreadLocali, Li);
		const char *Labeli = InputSeqs.GetLabel(ThreadLocali);

		uint Lj;
		const byte *Seqj = InputSeqs.GetByteSeq(ThreadLocalj, Lj);
		const char *Labelj = InputSeqs.GetLabel(ThreadLocalj);

		PathInfo *PI = ObjMgr::GetPathInfo();
		XDPMem &Mem = GetDPMem();
		ViterbiFastMem(Mem, Seqi, Li, Seqj, Lj, *PI);

		string RowX, RowY;
		MakeAlnRows(Seqi, Li, Seqj, Lj, *PI, RowX, RowY);
		const uint ColCount = PI->GetColCount();
		asserta(SIZE(RowX) == ColCount);
		asserta(SIZE(RowY) == ColCount);
		ObjMgr::Down(PI);

		double dij = GetProtDist(RowX.c_str(), RowY.c_str(), ColCount);
		if (fOut != 0)
			{
			LOCK();
			fprintf(fOut, "%s\t%s\t%.4g\n", Labeli, Labelj, dij);
			UNLOCK();
			}
		}
	}
