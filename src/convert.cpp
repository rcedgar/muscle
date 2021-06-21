#include "muscle.h"
#include "seqvect.h"

Sequence *SeqToSequence(const Seq &inseq, uint Index)
	{
	const char *Label = inseq.GetName();
	Sequence *outseq = new Sequence;
	asserta(outseq != 0);
	const uint L = SIZE(inseq);
	//outseq->length = (int) L;
	vector<char> &data = *new vector<char>;
	data.resize(L+1);
	data[0] = '@';
	for (uint i = 0; i < L; ++i)
		data[i+1] = inseq[i];
	outseq->data = &data;
	outseq->label = mystrsave(Label);
	return outseq;
	}

void GetLabelToIndex(const SeqVect &SV, map<string, uint> &LabelToIndex)
	{
	LabelToIndex.clear();
	const uint SeqCount = SV.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = string(SV.GetSeqName(SeqIndex));
		if (LabelToIndex.find(Label) != LabelToIndex.end())
			Die("Duplicate label >%s", Label.c_str());
		LabelToIndex[Label] = SeqIndex;
		}
	}

void GetSeqsByLabels(const SeqVect &SV, const vector<string> &Labels,
  const map<string, uint> &LabelToIndex, MultiSequence &Seqs)
	{
	Seqs.sequences = new vector<Sequence*>;
	Seqs.sequences->clear();
	const uint n = SIZE(Labels);
	for (uint i = 0; i < n; ++i)
		{
		const string &Label = Labels[i];
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			Die("Label not found >%s", Label.c_str());
		uint SeqIndex = p->second;
		const Seq &MySeq = SV.GetSeq(SeqIndex);
		Sequence *seq = SeqToSequence(MySeq, SeqIndex);
		Seqs.AddSequence(seq);
		}
	}
