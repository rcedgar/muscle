#include "muscle.h"
#include "mpcflat.h"

float CalcAlnFlat(const float *Post, uint LX, uint LY,
  float *DPRows, char *TB, string &Path);

MultiSequence *SqueezeGappyCols(const MultiSequence &Aln);

MultiSequence *MPCFlat::AlignAlns(const MultiSequence *ptrMSA1,
  const MultiSequence *ptrMSA2, float *ptrScore)
	{
	const uint MAX_COL_COUNT = optd(maxcols, 5000);
	const uint SeqCount1 = ptrMSA1->GetSeqCount();
	const uint SeqCount2 = ptrMSA2->GetSeqCount();

	uint ColCount1 = ptrMSA1->GetColCount();
	uint ColCount2 = ptrMSA2->GetColCount();

	MultiSequence *ptrSqueezedMSA1 = 0;
	MultiSequence *ptrSqueezedMSA2 = 0;
	if (opt(squeeze) && ColCount1 > MAX_COL_COUNT)
		{
		ptrMSA1 = ptrSqueezedMSA1 = SqueezeGappyCols(*ptrMSA1);
		ColCount1 = ptrMSA1->GetColCount();
		}

	if (opt(squeeze) && ColCount2 > MAX_COL_COUNT)
		{
		ptrMSA2 = ptrSqueezedMSA2 = SqueezeGappyCols(*ptrMSA2);
		ColCount2 = ptrMSA2->GetColCount();
		}

// INT_MAX=2147483647; 5 N^2 == 2147483647; N = sqrt(2147483647/5) = 20724
	if (double(ColCount1)*double(ColCount2)*5 + 100 > double(INT_MAX))
		Die("Join Cols1=%u, Cols2=%u overflow 64-bit HMM buffers", ColCount1, ColCount2);

	float *Post = AllocPost(ColCount1, ColCount2);
	BuildPost(*ptrMSA1, *ptrMSA2, Post);

	float *DPRows = AllocDPRows(ColCount1, ColCount2);
	char *TB = AllocTB(ColCount1, ColCount2);

	string Path;
	float Score = CalcAlnFlat(Post, ColCount1, ColCount2, DPRows, TB, Path);
	if (ptrScore != 0)
		*ptrScore = Score;
	myfree(Post);
	myfree(DPRows);
	myfree(TB);

	MultiSequence *result = new MultiSequence();
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount1; ++SeqIndex1)
		{
		const Sequence *InputRow = ptrMSA1->GetSequence(SeqIndex1);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'X');
		result->AddSequence(AlignedRow, true);
		}

	for (uint SeqIndex2 = 0; SeqIndex2 < SeqCount2; ++SeqIndex2)
		{
		const Sequence *InputRow = ptrMSA2->GetSequence(SeqIndex2);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'Y');
		result->AddSequence(AlignedRow, true);
		}

	if (ptrSqueezedMSA1 != 0) delete ptrSqueezedMSA1;
	if (ptrSqueezedMSA2 != 0) delete ptrSqueezedMSA2;

	return result;
	}
