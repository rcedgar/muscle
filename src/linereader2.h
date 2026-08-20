#pragma once
#include <mutex>

const unsigned LR_BUFF2 = (32*1024*1024);

// Returns strings instead of t_LineBuff
class LineReader2
	{
	mutex m_Lock;

public:
	string m_FileName;
	FILE *m_f;
	unsigned m_LineNr;
	char *m_Buffer;
	uint32 m_BufferBytes;
	uint32 m_BufferOffset;
	uint64 m_FileSize;
	bool m_EOF;

public:
	LineReader2();
	virtual ~LineReader2();

public:
	void Open(const string &FileName);
	void Close();
	void Rewind();
	unsigned GetPctDoneX10();
	double GetPctDone();

// Caller must own memory for Line because
// LineReader object may be used by multiple threads.
	bool ReadLine(string &Line);
	uint64 GetPos() const;

protected:
	void FillBuff();
	void Clear();
	};
