#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

const unsigned FASTA_BLOCK = 60;

void MSA::FromFASTAFile(TextFile &File)
	{
	Clear();

	FILE *f = File.GetStdioFile();
	
	unsigned uSeqCount = 0;
	unsigned uColCount = uInsane;
	for (;;)
		{
		char *Label;
		unsigned uSeqLength;
		char *SeqData = GetFastaSeq(f, &uSeqLength, &Label, false);
		if (0 == SeqData)
			break;
		AppendSeq(SeqData, uSeqLength, Label);
		}
	}

void MSA::FromFASTAFile_PreserveCase(const string &FileName)
	{
	extern bool g_FASTA_Upper;
	bool SaveUpper = g_FASTA_Upper;
	g_FASTA_Upper = false;
	FromFASTAFile(FileName);
	g_FASTA_Upper = true;
	}

void MSA::FromFASTAFile(const string &FileName)
	{
	TextFile TF(FileName);
	FromFASTAFile(TF);
	TF.Close();
	}

void MSA::ToFASTAFile(const string &FileName) const
	{
	TextFile TF(FileName, true);
	ToFASTAFile(TF);
	TF.Close();
	}

void MSA::ToFASTAFile(TextFile &File) const
	{
	const unsigned uColCount = GetColCount();
	assert(uColCount > 0);
	const unsigned uLinesPerSeq = (GetColCount() - 1)/FASTA_BLOCK + 1;
	const unsigned uSeqCount = GetSeqCount();

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		File.PutString(">");
		File.PutString(GetSeqName(uSeqIndex));
		File.PutString("\n");

		unsigned n = 0;
		for (unsigned uLine = 0; uLine < uLinesPerSeq; ++uLine)
			{
			unsigned uLetters = uColCount - uLine*FASTA_BLOCK;
			if (uLetters > FASTA_BLOCK)
				uLetters = FASTA_BLOCK;
			for (unsigned i = 0; i < uLetters; ++i)
				{
				char c = GetChar(uSeqIndex, n);
				File.PutChar(c);
				++n;
				}
			File.PutChar('\n');
			}
		}
	}
