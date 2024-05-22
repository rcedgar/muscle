#include "muscle.h"

void AlignTwoMSAsGivenPath(const MultiSequence &msaA,
  const MultiSequence &msaB, const string &Path,
  MultiSequence &msa2)
	{
	msa2.Clear();
	const uint SeqCountA = msaA.GetSeqCount();
	const uint SeqCountB = msaB.GetSeqCount();
	for (uint SeqIndexA = 0; SeqIndexA < SeqCountA; ++SeqIndexA)
		{
		const Sequence *InputRow = msaA.GetSequence(SeqIndexA);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'D');
		msa2.AddSequence(AlignedRow, true);
		}

	for (uint SeqIndexB = 0; SeqIndexB < SeqCountB; ++SeqIndexB)
		{
		const Sequence *InputRow = msaB.GetSequence(SeqIndexB);
		Sequence *AlignedRow = InputRow->AddGapsPath(Path, 'I');
		msa2.AddSequence(AlignedRow, true);
		}
	}
