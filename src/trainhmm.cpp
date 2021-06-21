#include "muscle.h"

static float FixP(float P)
	{
	if (P < 0)
		P = 0;
	if (P > 1)
		P = 1;
	return P;
	}

static void Train1(Sequence *seq1, Sequence *seq2,
  const vector<float> &forward, const vector<float> &backward,
  const float CurrentTransProb[HMMSTATE_COUNT][HMMSTATE_COUNT],
  const float CurrentMatchProb[256][256],
  const float CurrentInsProb[256],
  vector<float> &initDistribMat, vector<float> &gapOpen, vector<float> &gapExtend)
	{
	assert(seq1);
	assert(seq2);

	const int seq1Length = seq1->GetLength();
	const int seq2Length = seq2->GetLength();
	const char *iter1 = seq1->GetCharPtr1();
	const char *iter2 = seq2->GetCharPtr1();

	const float totalProb = PairHMM::ComputeTotalProbability(
	  seq1Length, seq2Length, forward, backward);

	vector<vector<float> > transCounts(HMMSTATE_COUNT, vector<float>(HMMSTATE_COUNT, LOG_ZERO));
	vector<float> initCounts(HMMSTATE_COUNT, LOG_ZERO);

	int ij = 0;
	int i1j = -seq2Length - 1;
	int ij1 = -1;
	int i1j1 = -seq2Length - 2;

	ij *= HMMSTATE_COUNT;
	i1j *= HMMSTATE_COUNT;
	ij1 *= HMMSTATE_COUNT;
	i1j1 *= HMMSTATE_COUNT;

	// compute initial distribution posteriors
	initCounts[0] = LOG_ADD(forward[0 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 1)] +
		backward[0 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 1)],
		forward[0 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)] +
		backward[0 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)]);

	for (int k = 0; k < InsertStateCount; k++)
		{
		initCounts[2 * k + 1] = LOG_ADD(forward[2 * k + 1 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 0)] +
			backward[2 * k + 1 + HMMSTATE_COUNT * (1 * (seq2Length + 1) + 0)],
			forward[2 * k + 1 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)] +
			backward[2 * k + 1 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)]);
		initCounts[2 * k + 2] = LOG_ADD(forward[2 * k + 2 + HMMSTATE_COUNT * (0 * (seq2Length + 1) + 1)] +
			backward[2 * k + 2 + HMMSTATE_COUNT * (0 * (seq2Length + 1) + 1)],
			forward[2 * k + 2 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)] +
			backward[2 * k + 2 + HMMSTATE_COUNT * ((seq1Length + 1) * (seq2Length + 1) - 1)]);
		}

// Sum to set expected counts in transCounts[][]
	for (int i = 0; i <= seq1Length; i++)
		{
		byte c1 = (i == 0) ? '~' : (byte)toupper(iter1[i]);

		for (int j = 0; j <= seq2Length; j++)
			{
			byte c2 = (j == 0) ? '~' : (byte)toupper(iter2[j]);

			if (i > 0 && j > 0)
				{
				for (int k = 0; k < HMMSTATE_COUNT; k++)
					{
					LOG_PLUS_EQUALS(transCounts[k][0],
						forward[k + i1j1] + CurrentTransProb[k][0] +
						CurrentMatchProb[c1][c2] + backward[0 + ij]);
					}
				}

			if (i > 0)
				{
				for (int k = 0; k < InsertStateCount; k++)
					{
					LOG_PLUS_EQUALS(transCounts[0][2 * k + 1],
						forward[0 + i1j] + CurrentTransProb[0][2 * k + 1] +
						CurrentInsProb[c1] + backward[2 * k + 1 + ij]);
				
					LOG_PLUS_EQUALS(transCounts[2 * k + 1][2 * k + 1],
						forward[2 * k + 1 + i1j] + CurrentTransProb[2 * k + 1][2 * k + 1] +
						CurrentInsProb[c1] + backward[2 * k + 1 + ij]);
					}
				}

			if (j > 0)
				{
				for (int k = 0; k < InsertStateCount; k++)
					{
					LOG_PLUS_EQUALS(transCounts[0][2 * k + 2],
						forward[0 + ij1] + CurrentTransProb[0][2 * k + 2] +
						CurrentInsProb[c2] + backward[2 * k + 2 + ij]);

					LOG_PLUS_EQUALS(transCounts[2 * k + 2][2 * k + 2],
						forward[2 * k + 2 + ij1] + CurrentTransProb[2 * k + 2][2 * k + 2] +
						CurrentInsProb[c2] + backward[2 * k + 2 + ij]);
					}
				}

			ij += HMMSTATE_COUNT;
			i1j += HMMSTATE_COUNT;
			ij1 += HMMSTATE_COUNT;
			i1j1 += HMMSTATE_COUNT;
			}
		}

// Divide by total prob
	for (int i = 0; i < HMMSTATE_COUNT; i++)
		{
		initCounts[i] -= totalProb;
		for (int j = 0; j < HMMSTATE_COUNT; j++)
			transCounts[i][j] -= totalProb;
		}

	float totalInitDistribCounts = 0;
	for (int i = 0; i < HMMSTATE_COUNT; i++)
		totalInitDistribCounts += exp(initCounts[i]);

	float P = FixP((float)exp(initCounts[0]) / totalInitDistribCounts);
	initDistribMat[0] += P;
	for (int k = 0; k < InsertStateCount; k++)
		{
		float val = (exp(initCounts[2 * k + 1]) + exp(initCounts[2 * k + 2])) / 2;
		float P = FixP(val / totalInitDistribCounts);
		initDistribMat[2 * k + 1] += P;
		initDistribMat[2 * k + 2] += P;
		}

	float inMatchStateCounts = 0;
	for (int i = 0; i < HMMSTATE_COUNT; i++)
		inMatchStateCounts += exp(transCounts[0][i]);

	for (int i = 0; i < InsertStateCount; i++)
		{
		float inGapStateCounts =
			exp(transCounts[2 * i + 1][0]) +
			exp(transCounts[2 * i + 1][2 * i + 1]) +
			exp(transCounts[2 * i + 2][0]) +
			exp(transCounts[2 * i + 2][2 * i + 2]);

		gapOpen[2 * i] += gapOpen[2 * i + 1] =
			(exp(transCounts[0][2 * i + 1]) +
			exp(transCounts[0][2 * i + 2])) /
			(2 * inMatchStateCounts);

		gapExtend[2 * i] += gapExtend[2 * i + 1] =
			(exp(transCounts[2 * i + 1][2 * i + 1]) +
			exp(transCounts[2 * i + 2][2 * i + 2])) /
			inGapStateCounts;
		}
	}

void cmd_trainhmm()
	{
	const string& InputFileName = opt(trainhmm);

	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName, true);

	InitProbcons();

	vector<float> gapOpen(2*InsertStateCount);
	vector<float> gapExtend(2*InsertStateCount);
	vector<float> initDistribMat(HMMSTATE_COUNT);

	const uint SeqCount = InputSeqs.GetSeqCount();
	uint PairIndex = 0;
	const uint PairCount = (SeqCount*(SeqCount + 1))/2;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		Sequence* Seq1 = InputSeqs.GetSequence(SeqIndex1);
		for (uint SeqIndex2 = 0; SeqIndex2 < SeqIndex1; ++SeqIndex2)
			{
			ProgressStep(PairIndex++, PairCount, "Aligning");

			Sequence* Seq2 = InputSeqs.GetSequence(SeqIndex2);

			vector<float>* forward =
				PairHMM::ComputeForwardMatrix(Seq1, Seq2); assert(forward);

			vector<float>* backward =
				PairHMM::ComputeBackwardMatrix(Seq1, Seq2); assert(backward);

			vector<float>* posterior =
				PairHMM::ComputePosteriorMatrix(Seq1, Seq2, *forward, *backward);

			Train1(Seq1, Seq2, *forward, *backward, 
			  PairHMM::m_TransScore, PairHMM::m_MatchScore, PairHMM::m_InsScore,
			  initDistribMat, gapOpen, gapExtend);

			delete forward;
			delete backward;
			delete posterior;
			}
		}

	asserta(PairCount > 0);
	for (uint i = 0; i < 2*InsertStateCount; ++i)
		{
		gapOpen[i] /= PairCount;
		gapExtend[i] /= PairCount;

		Log("%u open %.4g ext %.4g\n", i, gapOpen[i], gapExtend[i]);
		}

	for (uint i = 0; i < HMMSTATE_COUNT; ++i)
		{
		initDistribMat[i] /= PairCount;
		Log("%u init %.4g\n", i, initDistribMat[i]);
		}
	}
