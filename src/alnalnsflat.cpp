#include "muscle.h"
#include "mpcflat.h"

float CalcAlnFlat(const float *Post, uint LX, uint LY,
  float *DPRows, char *TB, string &Path);

MultiSequence *MPCFlat::AlignAlns(const MultiSequence &MSA1,
  const MultiSequence &MSA2, float *ptrScore)
	{
	const uint SeqCount1 = MSA1.GetSeqCount();
	const uint SeqCount2 = MSA2.GetSeqCount();

	const uint ColCount1 = MSA1.GetColCount();
	const uint ColCount2 = MSA2.GetColCount();

	float *Post = AllocPost(ColCount1, ColCount2);
	BuildPost(MSA1, MSA2, Post);

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
		const Sequence *InputRow = MSA1.GetSequence(SeqIndex1);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'X');
		result->AddSequence(AlignedRow, true);
		}

	for (uint SeqIndex2 = 0; SeqIndex2 < SeqCount2; ++SeqIndex2)
		{
		const Sequence *InputRow = MSA2.GetSequence(SeqIndex2);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'Y');
		result->AddSequence(AlignedRow, true);
		}

	return result;
	}
