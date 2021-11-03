#include "muscle.h"

const unsigned FASTA_BLOCK = 60;

void MSA::FromFASTAFile(TextFile &File)
	{
	Clear();

	FILE *f = File.GetStdioFile();
	
	unsigned uSeqCount = 0;
	unsigned uColCount = UINT_MAX;
	for (;;)
		{
		char *Label;
		unsigned uSeqLength;
		char *GetFastaSeq(FILE *f, unsigned *ptrSeqLength, char **ptrLabel, bool DeleteGaps);
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

void MSA::FromStrings(const vector<string> &Strings)
	{
	Clear();
	if (Strings.empty())
		Die("MSA::FromStrings, no data");

	vector<string> Labels;
	vector<string> Seqs;

	string CurrSeq;
	for (uint i = 0; i < SIZE(Strings); ++i)
		{
		const string &s = Strings[i];
		char s0 = s.c_str()[0];
		if (s0 == '>')
			{
			if (!Labels.empty())
				Seqs.push_back(CurrSeq);
			Labels.push_back(s.substr(1));
			CurrSeq.clear();
			}
		else
			{
			for (uint i = 0; i < SIZE(s); ++i)
				{
				char c = s[i];
				if (!isspace(c))
					CurrSeq.push_back(c);
				}
			}
		}
	Seqs.push_back(CurrSeq);
	FromStrings2(Labels, Seqs);
	}

void MSA::FromStrings2(const vector<string> &Labels, vector<string> &Seqs)
	{
	const uint SeqCount = SIZE(Labels);
	if (SIZE(Seqs) != SeqCount)
		Die("Invalid FASTA, %u labels %u seqs", SIZE(Labels), SIZE(Seqs));
	if (SeqCount == 0)
		Die("Empty FASTA");

	const uint ColCount = SIZE(Seqs[0]);
	SetSize(SeqCount, ColCount);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Label = Labels[SeqIndex].c_str();
		const string &Str = Seqs[SeqIndex];
		const uint n = SIZE(Str);
		if (n != ColCount)
			Die("MSA not aligned, seq lengths %u, %u", ColCount, n);
			  
		const char *S = Str.c_str();

		m_szNames[SeqIndex] = mystrsave(Label);
		m_szSeqs[SeqIndex] = mystrsave(S);
		}
	}

void MSA::FromFASTAFile(const string &FileName)
	{
	Clear();
	TextFile TF(FileName);
	FromFASTAFile(TF);
	TF.Close();
	}

void MSA::ToFASTAFile(const string &FileName) const
	{
	if (FileName.empty())
		return;
	TextFile TF(FileName, true);
	ToFASTAFile(TF);
	TF.Close();
	}

void MSA::ToFASTAFile(FILE *f) const
	{
	if (f == 0)
		return;
	for (uint SeqIndex = 0; SeqIndex < m_uSeqCount; ++SeqIndex)
		{
		const byte *S = (const byte *) m_szSeqs[SeqIndex];
		const char *Label = m_szNames[SeqIndex];
		SeqToFasta(f, S, m_uColCount, Label);
		}
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
