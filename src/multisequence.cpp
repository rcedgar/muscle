#include "muscle.h"
#include "alpha3.h"
#include "sort.h"

uint MultiSequence::m_NewCount;
uint MultiSequence::m_DeleteCount;

void MultiSequence::LoadMFA(const string& filename, bool stripGaps)
	{
	// try opening file
	FileBuffer infile(filename.c_str());

	if (infile.fail())
		Die("Cannot open %s, errno=%d %s",
		  filename.c_str(), errno, strerror(errno));

	// if successful, then load using other LoadMFA() routine
	LoadMFA(infile, stripGaps);

	infile.close();
	}

#if SEQ_TRACE
void MultiSequence::AssertSeqIds() const
	{
	const uint SeqCount = GetSeqCount();
	asserta(SIZE(m_Owners) == SIZE(m_Seqs));
	for (uint i = 0; i < SeqCount; ++i)
		m_Seqs[i]->AssertId();
	}
#endif

void MultiSequence::Clear()
	{
	const uint SeqCount = GetSeqCount();
	asserta(SIZE(m_Owners) == SIZE(m_Seqs));
	for (uint i = 0; i < SeqCount; ++i)
		{
		if (m_Owners[i])
			DeleteSequence(m_Seqs[i]);
		}
	m_Seqs.clear();
	m_Owners.clear();
	}

void MultiSequence::Copy(const MultiSequence &rhs)
	{
	Clear();
	const uint SeqCount = rhs.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *r = rhs.GetSequence(i)->Clone();
		AddSequence(r, true);
		}
	}

bool MultiSequence::ColIsAllGaps(uint Col) const
	{
	const uint SeqCount = GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = GetChar(SeqIndex, Col);
		if (!isgap(c))
			return false;
		}
	return true;
	}

void MultiSequence::LogMe() const
	{
	Log("\n");
	Log("MultiSequence::LogMe(%p), %u seqs\n", this, GetSeqCount());
	for (uint i = 0; i < GetSeqCount(); ++i)
		m_Seqs[i]->LogMe();
	}

uint MultiSequence::GetSeqLength(uint SeqIndex) const
	{
	const Sequence *Seq = GetSequence(SeqIndex);
	uint L = (uint) Seq->GetLength();
	return L;
	}

void MultiSequence::GetLengthOrder(vector<uint> &SeqIndexes) const
	{
	const uint SeqCount = GetSeqCount();
	vector<uint> Ls;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint L = (uint) GetSequence(SeqIndex)->GetLength();
		Ls.push_back(L);
		}
	SeqIndexes.resize(SeqCount);
	const uint *PtrLs = Ls.data();
	uint *PtrOrder = SeqIndexes.data();
	QuickSortOrderDesc<uint>(PtrLs, SeqCount, PtrOrder);
	}

void MultiSequence::LoadMFA(FileBuffer& infile, bool stripGaps)
	{
	if (infile.fail())
		Die("LoadMFA read error");

	set<string> Labels;
	unsigned DupeCount = 0;
	for (;;)
		{
		Sequence *seq = Sequence::NewSequence();
		bool Ok = seq->FromFileBuffer(infile, stripGaps);
		if (!Ok)
			{
			DeleteSequence(seq);
			break;
			}

		string Label = seq->m_Label;
		bool Dupe = false;
		for (uint i = 1; i < 100; ++i)
			{
			if (Labels.find(Label) == Labels.end())
				break;
			if (!m_DupeLabelsOk)
				{
				Dupe = true;
				Ps(Label, "%s dupelabel%u", seq->m_Label.c_str(), i);
				}
			}
		if (Dupe)
			{
			Log("Duplicate label >%s", seq->m_Label.c_str());
			seq->m_Label = Label;
			++DupeCount;
			}
		Labels.insert(Label);
		m_Seqs.push_back(seq);
		m_Owners.push_back(true);
		}
	if (DupeCount > 0 && !m_DupeLabelsOk)
		Warning("%u duplicate labels", DupeCount);
	}

uint MultiSequence::GetColCount() const
	{
	assert(IsAligned());
	const Sequence *Seq = GetSequence(0);
	uint L = Seq->GetLength();
	return L;
	}

bool MultiSequence::IsAligned() const
	{
	int N = GetNumSequences();
	if (N == 0)
		return false;
	uint ColCount0 = GetSequence(0)->GetLength();
	for (int i = 1; i < N; ++i)
		{
		int ColCount = GetSequence(i)->GetLength();
		if (ColCount != ColCount0)
			return false;
		}
	return true;
	}

uint MultiSequence::GetSeqIndex(const string &Label, bool FailOnError) const
	{
	int N = GetNumSequences();
	for (int i = 0; i < N; ++i)
		{
		const Sequence *Seq = GetSequence(i);
		if (Seq->GetLabel() == Label)
			return i;
		}
	if (FailOnError)
		Die("Label not found >%s", Label.c_str());
	return UINT_MAX;
	}

bool MultiSequence::GuessIsNucleo() const
	{
// If at least MIN_NUCLEO_PCT of the first CHAR_COUNT non-gap
// letters belong to the nucleotide alphabet, guess nucleo.
// Otherwise amino.
	const unsigned CHAR_COUNT = 100;
	const unsigned MIN_NUCLEO_PCT = 95;

	const unsigned SeqCount = GetSeqCount();
	uint NucleoCount = 0;
	for (uint i = 0; i < 100; ++i)
		{
		uint SeqIndex = randu32()%SeqCount;
		const Sequence &seq = *GetSequence(SeqIndex);
		const uint L = seq.GetLength();
		uint r = randu32();
		const uint Pos = r%L;
		byte c = (byte) GetChar(SeqIndex, Pos);
		uint Letter = g_CharToLetterNucleo[c];
		if (Letter < 4)
			++NucleoCount;
		}
	if (NucleoCount > 75)
		return true;
	return false;
	}

void MultiSequence::FromStrings(const vector<string> &Labels,
  const vector<string> &Seqs)
	{
	Clear();
	const uint N = SIZE(Seqs);
	asserta(SIZE(Labels) == N);
	for (uint i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		const string &Str = Seqs[i];
		Sequence *Seq = Sequence::NewSequence();
		Seq->FromString(Label, Str);
		AddSequence(Seq, true);
		}
	}

void MultiSequence::ToMSA(MSA &msa) const
	{
	uint SeqCount = GetSeqCount();
	uint ColCount = GetColCount();
	msa.SetSize(SeqCount, ColCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const Sequence *Seq = GetSequence(SeqIndex);
		const byte *ByteSeq = Seq->GetBytePtr();
		const char *Label = Seq->GetLabel().c_str();
		uint L = Seq->GetLength();

		char *CharSeq = myalloc(char, L+1);
		memcpy(CharSeq, ByteSeq, L);
		CharSeq[L] = 0;
		msa.m_szSeqs[SeqIndex] = CharSeq;
		msa.m_szNames[SeqIndex] = mystrsave(Label);
		}
	}

double MultiSequence::GetMeanSeqLength() const
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount == 0)
		return 0;
	double SumSeqLength = 0;
	for (uint i = 0; i < SeqCount; ++i)
		SumSeqLength += GetSequence(i)->GetLength();
	return SumSeqLength/SeqCount;
	}

uint MultiSequence::GetMaxSeqLength() const
	{
	const uint SeqCount = GetSeqCount();
	uint MaxSeqLength = 0;
	for (uint i = 0; i < SeqCount; ++i)
		MaxSeqLength = max(MaxSeqLength, GetSequence(i)->GetLength());
	return MaxSeqLength;
	}

uint MultiSequence::GetMinSeqLength() const
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount == 0)
		return 0;
	uint MinSeqLength = UINT_MAX;
	for (uint i = 0; i < SeqCount; ++i)
		MinSeqLength = min(MinSeqLength, GetSequence(i)->GetLength());
	return MinSeqLength;
	}
