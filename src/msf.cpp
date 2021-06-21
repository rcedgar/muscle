#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

const int MAX_NAME = 63;

const unsigned uCharsPerLine = 50;
const unsigned uCharsPerBlock = 10;

// Truncate at first white space or MAX_NAME, whichever comes
// first, then pad with blanks up to PadLength.
static const char *GetPaddedName(const char *Name, int PadLength)
	{
	static char PaddedName[MAX_NAME+1];
	memset(PaddedName, ' ', MAX_NAME);
	size_t n = strcspn(Name, " \t");
	memcpy(PaddedName, Name, n);
	PaddedName[PadLength] = 0;
	return PaddedName;
	}

static const char *strfind(const char *s, const char *t)
	{
	size_t n = strcspn(s, t);
	if (0 == n)
		return 0;
	return s + n;
	}

// GCG checksum code kindly provided by Eric Martel.
unsigned MSA::GetGCGCheckSum(unsigned uSeqIndex) const
	{
	unsigned CheckSum = 0;
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		unsigned c = (unsigned) GetChar(uSeqIndex, uColIndex);
		CheckSum += c*(uColIndex%57 + 1);
		CheckSum %= 10000;		
		}
	return CheckSum;
	}

static void MSFFixGaps(MSA &a)
	{
	const int SeqCount = a.GetSeqCount();
	const int ColCount = a.GetColCount();
	for (int SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		for (int ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			if (a.IsGap(SeqIndex, ColIndex))
				a.SetChar(SeqIndex, ColIndex, '.');
		}
	}

void MSA::ToMSFFile(TextFile &File, const char *ptrComment) const
	{
// Cast away const, yuck
	SetMSAWeightsMuscle((MSA &) *this);
	MSFFixGaps((MSA &) *this);

	File.PutString("PileUp\n");
	
	if (0 != ptrComment)
		File.PutFormat("Comment: %s\n", ptrComment);
	else
		File.PutString("\n");

	char seqtype = (g_Alpha == ALPHA_DNA || g_Alpha == ALPHA_RNA) ? 'N' : 'P';
	File.PutFormat("  MSF: %u  Type: %c  Check: 0000  ..\n\n",
	  GetColCount(), seqtype);

	int iLongestNameLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char *Name = GetSeqName(uSeqIndex);
		const char *PaddedName = GetPaddedName(Name, MAX_NAME);
		int iLength = (int) strcspn(PaddedName, " \t");
		if (iLength > iLongestNameLength)
			iLongestNameLength = iLength;
		}
		
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char *Name = GetSeqName(uSeqIndex);
		const char *PaddedName = GetPaddedName(Name, iLongestNameLength);
		File.PutFormat(" Name: %s", PaddedName);
		File.PutFormat("  Len: %u  Check: %5u  Weight: %g\n",
		  GetColCount(), GetGCGCheckSum(uSeqIndex), GetSeqWeight(uSeqIndex));
		}
	File.PutString("\n//\n");
	if (0 == GetColCount())
		return;

	unsigned uLineCount = (GetColCount() - 1)/uCharsPerLine + 1;
	for (unsigned uLineIndex = 0; uLineIndex < uLineCount; ++uLineIndex)
		{
		File.PutString("\n");
		unsigned uStartColIndex = uLineIndex*uCharsPerLine;
		unsigned uEndColIndex = uStartColIndex + uCharsPerLine - 1;
		if (uEndColIndex >= GetColCount())
			uEndColIndex = GetColCount() - 1;
		for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
			{
			const char *Name = GetSeqName(uSeqIndex);
			const char *PaddedName = GetPaddedName(Name, iLongestNameLength);
			File.PutFormat("%s   ", PaddedName);
			for (unsigned uColIndex = uStartColIndex; uColIndex <= uEndColIndex;
			  ++uColIndex)
				{
				if (0 == uColIndex%uCharsPerBlock)
					File.PutString(" ");
				char c = GetChar(uSeqIndex, uColIndex);
				File.PutFormat("%c", c);
				}
			File.PutString("\n");
			}
		}
	}
