#include "muscle.h"
#include "alpha3.h"
#include "sort.h"

MultiSequence::~MultiSequence()
	{
	// if sequences allocated
	if (sequences)
		{
		// free all sequences
		for (vector<Sequence*>::iterator iter = sequences->begin(); iter != sequences->end(); ++iter)
			{
			assert(*iter);
			delete* iter;
			*iter = NULL;
			}

		// free sequence vector
		delete sequences;
		sequences = NULL;
		}
	}

// Call before d'tor when sequences were not allocated by this object
void MultiSequence::FreeNonOwner()
	{
	if (sequences == 0)
		return;
	sequences->clear();
	delete sequences;
	sequences = 0;
	}

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

void MultiSequence::Copy(const MultiSequence &rhs)
	{
	asserta(sequences == 0);
	const uint SeqCount = rhs.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		{
		Sequence *newseq = new Sequence;
		asserta(newseq);
		const Sequence *rhsseq = rhs.GetSequence(i);
		newseq->Copy(*rhsseq);
		AddSequence(newseq);
		}
	}

uint MultiSequence::GetGSI(uint SeqIndex) const
	{
	const Sequence *Seq = GetSequence(SeqIndex);
	uint L = (uint) Seq->GetGSI();
	return L;
	}

uint MultiSequence::GetLength(uint SeqIndex) const
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

void MultiSequence::LogGSIs(const char *Msg) const
	{
	const MultiSequence &GlobalMS = GetGlobalInputMS();
	Log("\n");
	Log("LogGSIs()");
	if (Msg != 0)
		Log("  %s  ", Msg);
	const uint SeqCount = GetSeqCount();
	Log("%u seqs\n", SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		Log("[%5u]", i);
		const Sequence *Seq = GetSequence(i);
		uint GSI = (uint) Seq->GetGSI();
		const string Label = string(Seq->GetLabel());
		const string GlobalLabel = string(GlobalMS.GetLabel(GSI));
		Log("  %5d", GSI);
		Log("  >%s", Label.c_str());
		Log("  (%s)", GlobalLabel.c_str());
		if (Label != GlobalLabel)
			Log("  <<< ERROR");
		Log("\n");
		}
	}

void MultiSequence::LoadMFA(FileBuffer& infile, bool stripGaps)
	{
	if (infile.fail())
		Die("LoadMFA read error");

	for (;;)
		{
		Sequence* seq = new Sequence;
		bool Ok = seq->FromFileBuffer(infile, stripGaps);
		if (!Ok)
			{
			delete seq;
			break;
			}

		if (sequences == 0)
			sequences = new vector<Sequence *>;
		sequences->push_back(seq);
		}
	if (!sequences)
		Die("No sequences");
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
