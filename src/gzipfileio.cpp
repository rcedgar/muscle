#include "myutils.h"
#undef brk
#if 1 // _MSC_VER
#include "zlib.h"
#else
#include <zlib.h>
#endif

FILE *OpenGzipFile(const string &FileName)
	{
	gzFile f = gzopen(FileName.c_str(), "rb");
	if (f == 0)
		Die("Error opening gzip file %s", FileName.c_str());
	return (FILE *) f;
	}

uint32 ReadGzipFile(FILE *f, void *Buff, uint32 MaxBytes)
	{
	int n = gzread(gzFile(f), Buff, MaxBytes);
	if (n < 0)
		Die("Error reading gzip file");
	return unsigned(n);
	}

void RewindGzipFile(FILE *f)
	{
	int rc = gzrewind(gzFile(f));
	if (rc < 0)
		Die("gzrewind=%d", rc);
	}

uint64 GetGzipFileSize_NoFail(FILE *f)
	{
	long CurrPos = gzseek(gzFile(f), 0, 1);
	if (CurrPos < 0)
		return UINT64_MAX;
	long FileSize = gzseek(gzFile(f), 0, 2);
	uint64 Size64 = (FileSize >= 0 ? uint64(FileSize) : UINT64_MAX);
	gzseek(gzFile(f), CurrPos, 0);
	return Size64;
	}

uint64 GetGzipFilePos(FILE *f)
	{
	uint64 Pos = gzseek(gzFile(f), 0, 1);
	return Pos;
	}

void CloseGzipFile(FILE *f)
	{
	if (f == 0)
		return;
	gzclose_r(gzFile(f));
	}

void ReadLinesFromGzipFile(const string &FileName, vector<string> &Lines)
	{
	Lines.clear();
	FILE *f = OpenGzipFile(FileName);
	const uint BUFFER_SIZE = 1024*1024;
	byte *Data = myalloc(byte, BUFFER_SIZE);
	string Line;
	for (;;)
		{
		uint BytesRead = ReadGzipFile(f, Data, BUFFER_SIZE);
		if (BytesRead == 0)
			{
			if (!Line.empty())
				Lines.push_back(Line);
			break;
			}
		for (uint i = 0; i < BytesRead; ++i)
			{
			char c = Data[i];
			if (c == '\r')
				continue;
			else if (c == '\n')
				{
				Lines.push_back(Line);
				Line.clear();
				}
			else
				Line.push_back(c);
			}
		}
	myfree(Data);
	CloseGzipFile(f);
	}

void cmd_gunzip()
	{
	const string &InputFileName = g_Arg1;
	const string &OutputFileName = opt(output);

	FILE *fIn = OpenGzipFile(InputFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	const unsigned M = 1024*1024;
	byte *Buffer = myalloc(byte, M);

	for (;;)
		{
		uint32 n = ReadGzipFile(fIn, Buffer, M);
		if (n == 0)
			break;
		WriteStdioFile(fOut, Buffer, n);
		}

	CloseGzipFile(fIn);
	CloseStdioFile(fOut);
	}

void cmd_gunzip_lines()
	{
	const string &InputFileName = g_Arg1;

	vector<string> Lines;
	ReadLinesFromGzipFile(InputFileName, Lines);

	if (!optset_output)
		return;

	FILE *fOut = CreateStdioFile(opt(output));
	for (uint i = 0; i < SIZE(Lines); ++i)
		{
		fputs(Lines[i].c_str(), fOut);
		fputc('\n', fOut);
		}
	CloseStdioFile(fOut);
	}
