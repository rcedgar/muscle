#include "muscle.h"
#include "textfile.h"
#include <errno.h>

TextFile::TextFile(const char szFileName[], bool bWrite)
	{
	FILE *ptrFile = 0;
	if (bWrite)
		{
		if (0 == strcmp(szFileName, "-"))
			ptrFile = stdout;
		else
			ptrFile = fopen(szFileName, "wb");
		}
	else
		{
		if (0 == strcmp(szFileName, "-"))
			ptrFile = stdin;
		else
			ptrFile = fopen(szFileName, "rb");
		}
	if (0 == ptrFile)
		Die("Cannot open '%s' errno=%d\n", szFileName, errno);
	Init(ptrFile, szFileName);
	}

TextFile::TextFile(const string &FileName, bool bWrite)
	{
	const char *szFileName = FileName.c_str();
	FILE *ptrFile = 0;
	if (bWrite)
		{
		if (0 == strcmp(szFileName, "-"))
			ptrFile = stdout;
		else
			ptrFile = fopen(szFileName, "wb");
		}
	else
		{
		if (0 == strcmp(szFileName, "-"))
			ptrFile = stdin;
		else
			ptrFile = fopen(szFileName, "rb");
		}
	if (0 == ptrFile)
		Die("Cannot open '%s' errno=%d\n", szFileName, errno);
	Init(ptrFile, szFileName);
	}

void TextFile::Init(FILE *ptrFile, const char *ptrFileName)
	{
	m_ptrFile = ptrFile;
	m_ptrName = strdup(ptrFileName);
	m_uLineNr = 1;
	m_uColNr = 0;
	m_bLastCharWasEOL = true;
	m_cPushedBack = -1;
#if	DEBUG
	setbuf(m_ptrFile, 0);
#endif
	}

TextFile::TextFile(FILE *ptrFile, const char *ptrFileName)
	{
	Init(ptrFile, "-");
	}

TextFile::~TextFile()
	{
	if (m_ptrFile &&
	  m_ptrFile != stdin && m_ptrFile != stdout && m_ptrFile != stderr)
		fclose(m_ptrFile);
	free(m_ptrName);
	}

// Get line from file.
// Return true if end-of-file, quit if line too long.
bool TextFile::GetLine(char szLine[], unsigned uBytes)
	{
	if (0 == uBytes)
		Die("TextFile::GetLine, buffer zero size");

	
	int FillVal = 0; // suppress warning from gcc that I don't understand
	memset(szLine, FillVal, (size_t) uBytes);

	unsigned uBytesCopied = 0;

// Loop until end of line or end of file.
	for (;;)
		{
		char c;
		bool bEof = GetChar(c);
		if (bEof)
			return true;
		if ('\r' == c)
			continue;
		if ('\n' == c)
			return false;
		if (uBytesCopied < uBytes - 1)
			szLine[uBytesCopied++] = (char) c;
		else
			Die("TextFile::GetLine: input buffer too small, line %u",
			  m_uLineNr);
		}
	}

// As GetLine, but trim leading and trailing blanks; skip empty lines
bool TextFile::GetTrimLine(char szLine[], unsigned uBytes)
	{
	Die("GetTrimLine");
	return false;
	}

void TextFile::Rewind()
	{
	fseek(m_ptrFile, 0, SEEK_SET);
	m_uLineNr = 1;
	m_bLastCharWasEOL = true;
	}

void TextFile::PutChar(char c)
	{
	int i = fputc(c, m_ptrFile);
	assert(i == c);
	if ('\n' == c)
		{
		++m_uLineNr;
		m_uColNr = 1;
		}
	else
		++m_uColNr;
	}

void TextFile::PutString(const char szLine[])
	{
	int iError = fputs(szLine, m_ptrFile);
	assert(iError >= 0);
	}

void TextFile::PutFormat(const char szFormat[], ...)
	{
	char szStr[4096];
	va_list ArgList;
	va_start(ArgList, szFormat);
	vsprintf(szStr, szFormat, ArgList);
	PutString(szStr);
	}

void TextFile::GetLineX(char szLine[], unsigned uBytes)
	{
	if (uBytes == 0)
		Die("GetLineX");
	bool bEof = GetLine(szLine, uBytes);
	if (bEof)
		Die("end-of-file in GetLineX");
	}

bool TextFile::GetToken(char szToken[], unsigned uBytes, const char szCharTokens[])
	{
// Skip leading white space
	char c;
	for (;;)
		{
		bool bEof = GetChar(c);
		if (bEof)
			return true;
		if (!isspace(c))
			break;
		}

// Check for special case single-character tokens
	if (0 != strchr(szCharTokens, c))
		{
		assert(uBytes >= 2);
		szToken[0] = c;
		szToken[1] = 0;
		return false;
		}

// Loop until token terminated by white space, EOF or special
	unsigned uBytesCopied = 0;
	for (;;)
		{
		if (uBytesCopied < uBytes - 1)
			szToken[uBytesCopied++] = c;
		else
			Die("TextFile::GetToken: input buffer too small, line %u",
			  m_uLineNr);
		bool bEof = GetChar(c);
		if (bEof)
			{
			szToken[uBytesCopied] = 0;
			return true;
			}
	// Check for special case single-character tokens
		if (0 != strchr(szCharTokens, c))
			{
			PushBack(c);
			assert(uBytesCopied > 0 && uBytesCopied < uBytes);
			szToken[uBytesCopied] = 0;
			return false;
			}
		if (isspace(c))
			{
			assert(uBytesCopied > 0 && uBytesCopied < uBytes);
			szToken[uBytesCopied] = 0;
			return false;
			}
		}
	}

void TextFile::GetTokenX(char szToken[], unsigned uBytes, const char szCharTokens[])
	{
	bool bEof = GetToken(szToken, uBytes, szCharTokens);
	if (bEof)
		Die("End-of-file in GetTokenX");
	}

void TextFile::Skip()
	{
	for (;;)
		{
		char c;
		bool bEof = GetChar(c);
		if (bEof || '\n' == c)
			return;
		assert(isspace(c));
		}
	}

#ifdef _WIN32

TEXTFILEPOS TextFile::GetPos()
	{
	fpos_t p;
	int i = fgetpos(m_ptrFile, &p);
	assert(0 == i);
	assert(p >= 0);
	TEXTFILEPOS Pos;
	Pos.uOffset = (unsigned) p;
	Pos.uLineNr = m_uLineNr;
	Pos.uColNr = m_uColNr;
	return Pos;
	}

void TextFile::SetPos(TEXTFILEPOS Pos)
	{
	fpos_t p = (fpos_t) Pos.uOffset;
	int i = fsetpos(m_ptrFile, &p);
	assert(0 == i);
	m_uLineNr = Pos.uLineNr;
	m_uColNr = Pos.uColNr;
	}

#else

TEXTFILEPOS TextFile::GetPos()
	{
	TEXTFILEPOS Pos;
	Pos.uOffset = ftell(m_ptrFile);
	Pos.uLineNr = m_uLineNr;
	Pos.uColNr = m_uColNr;
	return Pos;
	}

void TextFile::SetPos(TEXTFILEPOS Pos)
	{
	fseek(m_ptrFile, Pos.uOffset, SEEK_SET);
	m_uLineNr = Pos.uLineNr;
	m_uColNr = Pos.uColNr;
	}

#endif

bool TextFile::GetChar(char &c)
	{
	if (-1 != m_cPushedBack)
		{
		c = (char) m_cPushedBack;
		m_cPushedBack = -1;
		return false;
		}

	int ic = fgetc(m_ptrFile);
	if (ic < 0)
		{
		if (feof(m_ptrFile))
			{
		// Hack to fix up a non-empty text file that is missing
		// and end-of-line character in the last line.
			if (!m_bLastCharWasEOL && m_uLineNr > 0)
				{
				c = '\n';
				m_bLastCharWasEOL = true;
				return false;
				}
			return true;
			}
		Die("TextFile::GetChar, error %s", strerror(errno));
		}
	c = (char) ic;
	if ('\n' == c)
		{
		m_bLastCharWasEOL = true;
		++m_uLineNr;
		m_uColNr = 1;
		}
	else
		{
		m_bLastCharWasEOL = false;
		++m_uColNr;
		}
	return false;
	}

void TextFile::GetCharX(char &c)
	{
	bool bEof = GetChar(c);
	if (bEof)
		Die("End-of-file in GetCharX");
	}

void TextFile::GetNonblankChar(char &c)
	{
	do
		{
		bool bEof = GetChar(c);
		if (bEof)
			Die("End-of-file in GetCharX");
		}
	while (isspace(c));
	}

void TextFile::SkipLine()
	{
	if (m_bLastCharWasEOL)
		return;
	for (;;)
		{
		char c;
		bool bEof = GetChar(c);
		if (bEof)
			Die("End-of-file in SkipLine");
		if ('\n' == c)
			break;
		}
	}

void TextFile::SkipWhite()
	{
	bool bEof = SkipWhiteX();
	if (bEof)
		Die("End-of-file skipping white space");
	}

bool TextFile::SkipWhiteX()
	{
	for (;;)
		{
		char c;
		bool bEof = GetChar(c);
		if (bEof)
			return true;
		if (!isspace(c))
			{
			PushBack(c);
			break;
			}
		}
	return false;
	}
