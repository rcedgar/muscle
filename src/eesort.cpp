#include "muscle.h"
#include "sort.h"
#include "locallock.h"

void cmd_eesort()
	{
	const string &QueryFileName = opt(eesort);
	const string &DBFileName = opt(db);
	const string &OutputFileName = opt(output);

	FILE *fTsv = CreateStdioFile(opt(tsvout));
	FILE *fFa = CreateStdioFile(OutputFileName);

	MultiSequence Query;
	MultiSequence DB;

	Query.FromFASTA(QueryFileName, true);

	Progress("Reading %s ...", DBFileName.c_str());
	DB.FromFASTA(DBFileName, true);
	Progress("done\n");

	bool IsNucleo = DB.GuessIsNucleo();
	if (IsNucleo)
		SetAlpha(ALPHA_Nucleo);
	else
		SetAlpha(ALPHA_Amino);
	InitProbcons();

	const uint QuerySeqCount = Query.GetSeqCount();
	const uint DBSeqCount = DB.GetSeqCount();

	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairCounter = 0;
	vector<double> EAs(DBSeqCount, DBL_MAX);

#pragma omp parallel for num_threads(ThreadCount)
	for (int iDBSeqIndex = 0; iDBSeqIndex < (int) DBSeqCount; ++iDBSeqIndex)
		{
		Lock();
		ProgressStep(PairCounter++, DBSeqCount, "Calculating");
		Unlock();
		uint DBSeqIndex = uint(iDBSeqIndex);

		const Sequence *DBSeq = DB.GetSequence(DBSeqIndex);
		const string &DBLabel = DBSeq->m_Label;

		for (uint QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
			{
			const Sequence *QSeq = Query.GetSequence(QuerySeqIndex);
			const string &QLabel = QSeq->m_Label;

			string Path;
			double EA = AlignPairFlat(QLabel, DBLabel, Path);
			if (QuerySeqIndex == 0)
				{
				Lock();
				EAs[DBSeqIndex] = EA;
				Unlock();
				}
			}
		}

	vector<uint> Order(DBSeqCount);
	QuickSortOrderDesc(EAs.data(), DBSeqCount, Order.data());

	for (uint k = 0; k < DBSeqCount; ++k)
		{
		ProgressStep(k, DBSeqCount, "Writing %s", OutputFileName.c_str());

		uint DBSeqIndex = Order[k];
		const Sequence *DBSeq = DB.GetSequence(DBSeqIndex);

		double EA = EAs[DBSeqIndex];
		asserta(EA != DBL_MAX);
		Pf(fTsv, "%.3g	%s\n", EA, DBSeq->GetLabel().c_str());

		DBSeq->WriteMFA(fFa);
		}
	}
