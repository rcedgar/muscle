#include "muscle.h"

float PairHMM::AlignPair(Sequence *Seq1, Sequence *Seq2,
  vector<char> *Path)
	{
	InitProbcons();

	int L1 = Seq1->GetLength();
	int L2 = Seq2->GetLength();
	asserta(L1 > 0);
	asserta(L2 > 0);

	vector<float>* forward =
	  PairHMM::ComputeForwardMatrix(Seq1, Seq2);
	asserta(forward != 0);

	vector<float>* backward =
	  PairHMM::ComputeBackwardMatrix(Seq1, Seq2);
	asserta(backward != 0);

	vector<float>* posterior =
	  PairHMM::ComputePosteriorMatrix(Seq1, Seq2, *forward, *backward);
	asserta(posterior != 0);

	pair<vector<char>*, float> alignment =
	  PairHMM::ComputeAlignment(L1, L2, *posterior);

	if (Path != 0)
		*Path = *alignment.first;
	float Score = alignment.second;
	float EA = Score/min(L1, L2);

	delete forward;
	delete backward;
	delete posterior;

	return EA;
	}

float PairHMM::AlignPair_StrPath(Sequence *Seq1, Sequence *Seq2,
  string &StrPath)
	{
	InitProbcons();

	StrPath.clear();

	vector<char> VecPath;
	float EA = AlignPair(Seq1, Seq2, &VecPath);
	const uint ColCount = SIZE(VecPath);
	StrPath.reserve(ColCount);
	for (uint i = 0; i < ColCount; ++i)
		{
		char c = VecPath[i];
		StrPath += c;
		}
	return EA;
	}
