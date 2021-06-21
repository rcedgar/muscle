#ifndef	TextFile_h
#define TextFile_h

#include <stdio.h>

struct TEXTFILEPOS
	{
	unsigned uOffset;
	unsigned uLineNr;
	unsigned uColNr;
	};

const unsigned TextFileBufferSize = 256;

class TextFile
	{
private:
// no default c'tor, not implemented
	TextFile();

public:
	virtual ~TextFile();

	TextFile(const char szFileName[], bool bWrite = false);
	TextFile(const string &FileName, bool bWrite = false);
	TextFile(FILE *ptrFile, const char *ptrFileName = "-");
	void Close() { fclose(m_ptrFile); m_ptrFile = 0; }

	bool GetLine(char szLine[], unsigned uBytes);
	bool GetTrimLine(char szLine[], unsigned uBytes);
	void GetLineX(char szLine[], unsigned uBytes);

	bool GetToken(char szToken[], unsigned uBytes, const char szCharTokens[] = "{}");
	void GetTokenX(char szToken[], unsigned uBytes, const char szCharTokens[] = "{}");

	void Skip();
	void SkipLine();
	void SkipWhite();
	bool SkipWhiteX();
	void Rewind();
	TEXTFILEPOS GetPos();
	void SetPos(TEXTFILEPOS Pos);
	bool GetChar(char &c);
	void GetCharX(char &c);
	void GetNonblankChar(char &c);

	unsigned GetLineNr() { return m_uLineNr; }

	void PutString(const char szLine[]);
	void PutFormat(const char szFormat[], ...);
	void PutChar(char c);

	const char *GetFileName() { return m_ptrName; }

	void PushBack(int c) { m_cPushedBack = c; }

	FILE *GetStdioFile() const { return m_ptrFile; }

private:
	void Init(FILE *ptrFile, const char *ptrFileName);

private:
	FILE *m_ptrFile;
	unsigned m_uLineNr;
	unsigned m_uColNr;
	char *m_ptrName;
	bool m_bLastCharWasEOL;
	int m_cPushedBack;
	};

#endif // TextFile_h
