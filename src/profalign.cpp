#include "muscle.h"
#include "locallock.h"

MultiSequence *MPCFlat::ProfAlign(const MultiSequence &MSA1,
  const MultiSequence &MSA2)
	{
	uint SeqCount1 = MSA1.GetSeqCount();
	uint ColCount1 = MSA1.GetColCount();
	uint SeqCount2 = MSA2.GetSeqCount();
	uint ColCount2 = MSA2.GetColCount();

	uint SMI = 0;
	MultiSequence CombinedSeqs;
	for (uint SeqIndex = 0; SeqIndex < SeqCount1; ++SeqIndex)
		{
		const Sequence *seq = MSA1.GetSequence(SeqIndex);
		Sequence *seq2 = seq->CopyDeleteGaps();
		CombinedSeqs.AddSequence(seq2, false);

		Sequence *HackSeq = (Sequence *) seq;
		HackSeq->m_SMI = SMI++;
		}

	for (uint SeqIndex = 0; SeqIndex < SeqCount2; ++SeqIndex)
		{
		const Sequence *seq = MSA2.GetSequence(SeqIndex);
		Sequence *seq2 = seq->CopyDeleteGaps();
		CombinedSeqs.AddSequence(seq2, false);

		Sequence *HackSeq = (Sequence *) seq;
		HackSeq->m_SMI = SMI++;
		}

	InitSeqs(&CombinedSeqs);
	InitPairs();
	uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	AllocPairCount(PairCount);
	InitDistMx();
	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairCounter = 0;
	uint PairCount2 = SeqCount1*SeqCount2;
#pragma omp parallel for num_threads(ThreadCount)
	for (int SeqIndex1 = 0; SeqIndex1 < (int) SeqCount1; ++SeqIndex1)
		{
		for (int SeqIndex2 = 0; SeqIndex2 < (int) SeqCount2; ++SeqIndex2)
			{
			Lock();
			ProgressStep(PairCounter++, PairCount2, "Calc posteriors");
			Unlock();

			uint PairIndex = GetPairIndex(SeqIndex1, SeqCount1 + SeqIndex2);
			CalcPosterior(PairIndex);
			}
		}

	float Score = 0;
	MultiSequence *MSA = AlignAlns(MSA1, MSA2, &Score);
	return MSA;
	}

void cmd_profalign()
	{
	MultiSequence MSA1;
	MultiSequence MSA2;
	MSA1.LoadMFA(opt(profalign), false);
	MSA2.LoadMFA(opt(input2), false);
	bool IsNucleo = MSA1.GuessIsNucleo();
	bool IsNucleo2 = MSA2.GuessIsNucleo();
	asserta(IsNucleo2 == IsNucleo);
	if (IsNucleo)
		SetAlpha(ALPHA_Nucleo);
	else
		SetAlpha(ALPHA_Amino);

	uint PerturbSeed = 0;
	if (optset_perturb)
		PerturbSeed = opt(perturb);

	HMMParams HP;
	HP.FromDefaults(IsNucleo);
	if (PerturbSeed > 0)
		{
		ResetRand(PerturbSeed);
		HP.PerturbProbs(PerturbSeed);
		}
	HP.ToPairHMM();

	MPCFlat M;
	MultiSequence *MSA = M.ProfAlign(MSA1, MSA2);
	MSA->WriteMFA(opt(output));
	}
 