#include "muscle.h"
#include "seqvect.h"
#include "textfile.h"
#include "msa.h"

const size_t MAX_FASTA_LINE = 16000;

SeqVect::~SeqVect()
	{
	Clear();
	}

void SeqVect::Clear()
	{
	for (size_t n = 0; n < size(); ++n)
		delete (*this)[n];
	}

void SeqVect::ToFASTAFile(const string &FileName) const
	{
	TextFile TF(FileName);
	ToFASTAFile(TF);
	TF.Close();
	}

void SeqVect::ToFASTAFile(TextFile &File) const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->ToFASTAFile(File);
		}
	}

void SeqVect::FromFASTAFile(const string &FileName)
	{
	TextFile TF(FileName);
	FromFASTAFile(TF);
	TF.Close();
	}

void SeqVect::FromFASTAFile(TextFile &File)
	{
	Clear();

	FILE *f = File.GetStdioFile();
	for (;;)
		{
		char *Label;
		unsigned uLength;
		char *SeqData = GetFastaSeq(f, &uLength, &Label);
		if (0 == SeqData)
			return;
		Seq *ptrSeq = new Seq;

		for (unsigned i = 0; i < uLength; ++i)
			{
			char c = SeqData[i];
			ptrSeq->push_back(c);
			}

		ptrSeq->SetName(Label);
		push_back(ptrSeq);

		delete[] SeqData;
		delete[] Label;
		}
	}

void SeqVect::PadToMSA(MSA &msa)
	{
	unsigned uSeqCount = Length();
	if (0 == uSeqCount)
		{
		msa.Clear();
		return;
		}

	unsigned uLongestSeqLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		unsigned uColCount = ptrSeq->Length();
		if (uColCount > uLongestSeqLength)
			uLongestSeqLength = uColCount;
		}
	msa.SetSize(uSeqCount, uLongestSeqLength);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		msa.SetSeqName(uSeqIndex, ptrSeq->GetName());
		unsigned uColCount = ptrSeq->Length();
		unsigned uColIndex;
		for (uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			char c = ptrSeq->at(uColIndex);
			msa.SetChar(uSeqIndex, uColIndex, c);
			}
		while (uColIndex < uLongestSeqLength)
			msa.SetChar(uSeqIndex, uColIndex++, '.');
		}
	}

void SeqVect::Copy(const SeqVect &rhs)
	{
	clear();
	unsigned uSeqCount = rhs.Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = rhs.at(uSeqIndex);
		Seq *ptrSeqCopy = new Seq;
		ptrSeqCopy->Copy(*ptrSeq);
		push_back(ptrSeqCopy);
		}
	}

void SeqVect::StripGaps()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGaps();
		}
	}

void SeqVect::StripGapsAndWhitespace()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGapsAndWhitespace();
		}
	}

void SeqVect::ToUpper()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->ToUpper();
		}
	}

bool SeqVect::FindName(const char *ptrName, unsigned *ptruIndex) const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		if (0 == stricmp(ptrSeq->GetName(), ptrName))
			{
			*ptruIndex = uSeqIndex;
			return true;
			}
		}
	return false;
	}

void SeqVect::AppendString(const string &Label, const string &s)
	{
	Seq &Copy = *new Seq;
	for (uint i = 0; i < SIZE(s); ++i)
		Copy.AppendChar(s[i]);

	Copy.SetName(Label.c_str());
	push_back(&Copy);
	}

void SeqVect::AppendSeq(const Seq &s)
	{
	Seq *ptrSeqCopy = new Seq;
	ptrSeqCopy->Copy(s);
	push_back(ptrSeqCopy);
	}

void SeqVect::LogMe() const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->LogMe();
		}
	}

void SeqVect::GetLabel(unsigned uSeqIndex, string &Label) const
	{
	const char *Name = GetSeqName(uSeqIndex);
	Label = string(Name);
	}

const char *SeqVect::GetSeqName(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetName();
	}

unsigned SeqVect::GetSeqId(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetId();
	}

unsigned SeqVect::GetSeqIdFromName(const char *Name) const
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned i = 0; i < uSeqCount; ++i)
		{
		if (!strcmp(Name, GetSeqName(i)))
			return GetSeqId(i);
		}
	Quit("SeqVect::GetSeqIdFromName(%s): not found", Name);
	return 0;
	}

Seq &SeqVect::GetSeqById(unsigned uId)
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned i = 0; i < uSeqCount; ++i)
		{
		if (GetSeqId(i) == uId)
			return GetSeq(i);
		}
	Quit("SeqVect::GetSeqIdByUd(%d): not found", uId);
	return (Seq &) *((Seq *) 0);
	}

unsigned SeqVect::GetSeqLength(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->Length();
	}

Seq &SeqVect::GetSeq(unsigned uSeqIndex)
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

const Seq &SeqVect::GetSeq(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

void SeqVect::SetSeqId(unsigned uSeqIndex, unsigned uId)
	{
	assert(uSeqIndex < size());
	Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->SetId(uId);
	}

ALPHA SeqVect::GuessAlpha() const
	{
// If at least MIN_NUCLEO_PCT of the first CHAR_COUNT non-gap
// letters belong to the nucleotide alphabet, guess nucleo.
// Otherwise amino.
	const unsigned CHAR_COUNT = 100;
	const unsigned MIN_NUCLEO_PCT = 95;

	const unsigned uSeqCount = GetSeqCount();
	if (0 == uSeqCount)
		return ALPHA_Amino;

	unsigned uSeqIndex = 0;
	unsigned uPos = 0;
	unsigned uSeqLength = GetSeqLength(0);
	unsigned uDNACount = 0;
	unsigned uRNACount = 0;
	unsigned uTotal = 0;
	const Seq *ptrSeq = &GetSeq(0);
	for (;;)
		{
		while (uPos >= uSeqLength)
			{
			++uSeqIndex;
			if (uSeqIndex >= uSeqCount)
				break;
			ptrSeq = &GetSeq(uSeqIndex);
			uSeqLength = ptrSeq->Length();
			uPos = 0;
			}
		if (uSeqIndex >= uSeqCount)
			break;
		char c = ptrSeq->at(uPos++);
		if (IsGapChar(c))
			continue;
		if (IsDNA(c))
			++uDNACount;
		if (IsRNA(c))
			++uRNACount;
		++uTotal;
		if (uTotal >= CHAR_COUNT)
			break;
		}
	if (uTotal != 0 && ((uDNACount*100)/uTotal) >= MIN_NUCLEO_PCT)
		return ALPHA_DNA;
	if (uTotal != 0 && ((uRNACount*100)/uTotal) >= MIN_NUCLEO_PCT)
		return ALPHA_RNA;
	return ALPHA_Amino;
	}

void SeqVect::FixAlpha()
	{
	ClearInvalidLetterWarning();
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->FixAlpha();
		}
	ReportInvalidLetters();
	}
