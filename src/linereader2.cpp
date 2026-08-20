#include "myutils.h"
#include "linereader2.h"

LineReader2::LineReader2()
	{
	Clear();
	}

LineReader2::~LineReader2()
	{
	Close();
	}

void LineReader2::Open(const string &FileName)
	{
	m_FileName = FileName;
	m_f = OpenStdioFile(FileName);
	m_Buffer = myalloc(char, LR_BUFF2);
	m_BufferBytes = 0;
	m_BufferOffset = 0;
	m_LineNr = 0;
	m_EOF = false;
	m_FileSize = GetStdioFileSize64(m_f);
	}

double LineReader2::GetPctDone()
	{
	if (m_FileSize == UINT64_MAX || m_FileSize == 0)
		return -1;
	uint64 Pos = GetPos();
	double Pct = double(Pos)/double(m_FileSize);
	return Pct;
	}

unsigned LineReader2::GetPctDoneX10()
	{
	if (m_FileSize == UINT64_MAX || m_FileSize == 0)
		return 0;
	uint64 Pos = GetPos();
	double f = double(Pos)/double(m_FileSize);
	uint n = uint(f*1000);
	if (n >= 999)
		n = 998;
	return n;
	}

void LineReader2::Clear()
	{
	m_f = 0;
	m_Buffer = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_LineNr = 0;
	m_EOF = true;
	}

void LineReader2::Close()
	{
	if (m_f == 0)
		return;

	CloseStdioFile(m_f);
	m_f = 0;
	myfree(m_Buffer);
	Clear();
	}

void LineReader2::Rewind()
	{
	SetStdioFilePos(m_f, 0);
	}

uint64 LineReader2::GetPos() const
	{
	return GetStdioFilePos64(m_f);
	}

bool LineReader2::ReadLine(string &Line)
	{
	Line.clear();
	if (m_EOF)
		return false;

	Line.reserve(1024);
	unsigned Length = 0;
	for (;;)
		{
		if (m_BufferOffset >= m_BufferBytes)
			{
			FillBuff();
			if (m_EOF)
				{
				if (Length == 0)
					return false;
				++m_LineNr;
				return true;
				}
			}

		char c = m_Buffer[m_BufferOffset++];
		if (c == '\r')
			continue;
		if (c == '\n')
			{
			++m_LineNr;
			return true;
			}
		Line += c;
		}
	}

void LineReader2::FillBuff()
	{
	if (m_EOF)
		return;
	uint32 BytesToRead = LR_BUFF2;
	m_BufferOffset = 0;
	m_BufferBytes = ReadStdioFile_NoFail(m_f, m_Buffer, LR_BUFF2);
	if (m_BufferBytes == 0)
		m_EOF = true;
	}
