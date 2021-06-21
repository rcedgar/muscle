#include "muscle.h"
#include "msa.h"
#include "textfile.h"

const int BLOCKSIZE = 60;

static char FixChar(char c)
	{
	switch (c)
		{
	case '(':
	case ')':
	case '[':
	case ']':
	case ':':
	case ';':
	case ',':
		return '_';
		}
	if (!isprint(c))
		return '_';
	return c;
	}

static void FixName(char Name[])
	{
	while (char c = *Name)
		*Name++ = FixChar(c);
	}

void MSA::ToPhySequentialFile(TextFile &File) const
	{
	const unsigned SeqCount = GetSeqCount();
	const unsigned ColCount = GetColCount();

	File.PutFormat("%d %d\n", SeqCount, ColCount);

	if (0 == ColCount)
		return;

	for (unsigned Seq = 0; Seq < SeqCount; ++Seq)
		{
		char Name[11];
		const char *ptrName = GetSeqName(Seq);
		size_t n = strlen(ptrName);
		if (n > 10)
			n = 10;
		memcpy(Name, ptrName, n);
		Name[n] = 0;
		FixName(Name);
		File.PutFormat("%-10.10s", Name);

		int BlockIndex = 0;
		int Col = 0;
		for (;;)
			{
			const unsigned MaxCols = (BlockIndex == 0) ? (BLOCKSIZE - 10) : BLOCKSIZE;
			for (unsigned ColsThisBlock = 0; ColsThisBlock < MaxCols; ++ColsThisBlock)
				{
				if (Col == ColCount)
					break;
				if (ColsThisBlock%10 == 0 && (BlockIndex == 0 || ColsThisBlock > 0))
					File.PutChar(' ');
				char c = GetChar(Seq, Col);
				if (isalpha(c))
					c = toupper(c);
				File.PutChar(c);
				++Col;
				}
			File.PutChar('\n');
			if (Col == ColCount)
				break;
			++BlockIndex;
			}
		}
	}

void MSA::ToPhyInterleavedFile(TextFile &File) const
	{
	const unsigned SeqCount = GetSeqCount();
	const unsigned ColCount = GetColCount();

	File.PutFormat("%d %d\n", SeqCount, ColCount);

	if (0 == ColCount)
		return;

	int Col = 0;
	for (;;)
		{
		const unsigned ColBlockStart = Col;
		const unsigned MaxCols = (ColBlockStart == 0) ? (BLOCKSIZE - 10) : BLOCKSIZE;

		for (unsigned Seq = 0; Seq < SeqCount; ++Seq)
			{
			if (0 == ColBlockStart)
				{
				char Name[11];
				const char *ptrName = GetSeqName(Seq);
				size_t n = strlen(ptrName);
				if (n > 10)
					n = 10;
				memcpy(Name, ptrName, n);
				Name[n] = 0;
				FixName(Name);
				File.PutFormat("%-10.10s", Name);
				}

			Col = ColBlockStart;
			for (unsigned ColsThisBlock = 0; ColsThisBlock < MaxCols; ++ColsThisBlock)
				{
				if (Col == ColCount)
					break;
				if (ColsThisBlock%10 == 0 && (0 == ColBlockStart || ColsThisBlock > 0))
					File.PutChar(' ');
				char c = GetChar(Seq, Col);
				if (isalpha(c))
					c = toupper(c);
				File.PutChar(c);
				++Col;
				}
			File.PutChar('\n');
			}
		if (Col == ColCount)
			break;
		File.PutChar('\n');
		}
	}
