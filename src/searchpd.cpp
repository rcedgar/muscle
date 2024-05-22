#include "muscle.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "omplock.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI);
void MakeAlnRows(const byte *XSeq, uint LX,
  const byte *YSeq, uint LY, const PathInfo &PI,
  string &RowX, string &RowY);
double GetProtDist(const char *Q, const char *T, uint ColCount);
XDPMem &GetDPMem();

double AlignAndProtDist(const byte *Seqi, uint Li,
  const byte *Seqj, uint Lj)
	{

	PathInfo *PI = ObjMgr::GetPathInfo();
	XDPMem &Mem = GetDPMem();
	ViterbiFastMem(Mem, Seqi, Li, Seqj, Lj, *PI);

	string RowX, RowY;
	MakeAlnRows(Seqi, Li, Seqj, Lj, *PI, RowX, RowY);

	const uint ColCount = PI->GetColCount();
	asserta(SIZE(RowX) == ColCount);
	asserta(SIZE(RowY) == ColCount);
	double dij = ::GetProtDist(RowX.c_str(), RowY.c_str(), ColCount);
	ObjMgr::Down(PI);
	return dij;
	}

void cmd_searchpd()
	{
	const string &InputFileName = opt(searchpd);
	asserta(optset_maxpd);
	asserta(optset_db);
	const string &DBFileName = opt(db);
	double MaxPD = opt(maxpd);
	if (optset_output)
		Die("Use -tsvout not -output");
	FILE *fOut = CreateStdioFile(opt(tsvout));

	MultiSequence Query;
	Query.FromFASTA(InputFileName, true);

	MultiSequence DB;
	DB.FromFASTA(DBFileName, true);

	const uint QuerySeqCount = Query.GetSeqCount();
	const uint DBSeqCount = DB.GetSeqCount();

	bool IsNucleo = Query.GuessIsNucleo();
	SetAlphab(IsNucleo);

	const uint SeqCount = Query.GetSeqCount();
	ProgressLog("%u query seqs, maxpd %.2f\n", SeqCount, MaxPD);

	uint Counter = 0;
	uint ThreadCount = GetRequestedThreadCount();
#pragma omp parallel for num_threads(ThreadCount)
	for (int QuerySeqIndex = 0; QuerySeqIndex < (int) QuerySeqCount;
	  ++QuerySeqIndex)
		{
		LOCK();
		ProgressStep(Counter++, QuerySeqCount, "Searching");
		UNLOCK();
		const Sequence *Q = Query.GetSequence(QuerySeqIndex);
		uint LQ = Q->GetLength();
		const byte *SeqQ = Q->GetBytePtr();

		for (uint DBSeqIndex = 0; DBSeqIndex < DBSeqCount;
		  ++DBSeqIndex)
			{
			const Sequence *T = DB.GetSequence(DBSeqIndex);
			uint LT = T->GetLength();
			const byte *SeqT = T->GetBytePtr();

			double d = AlignAndProtDist(SeqQ, LQ, SeqT, LT);
			if (d <= MaxPD)
				{
				if (fOut != 0)
					{
					LOCK();
					const char *LabelQ = Q->GetLabel().c_str();
					const char *LabelT = T->GetLabel().c_str();
					fprintf(fOut, "%s\t%s\t%.3g\n",
					  LabelQ, LabelT, d);
					UNLOCK();
					}
				}
			}
		}
	}
