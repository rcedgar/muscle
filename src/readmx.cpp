#include "muscle.h"
#include "textfile.h"

#define TRACE	0

const int MAX_LINE = 4096;
const int MAX_HEADINGS = 32;
static char Heading[MAX_HEADINGS];
static unsigned HeadingCount = 0;
static float Mx[32][32];

static void LogMx()
	{
	Log("Matrix\n");
	Log("     ");
	for (int i = 0; i < 20; ++i)
		Log("    %c", LetterToChar(i));
	Log("\n");

	for (int i = 0; i < 20; ++i)
		{
		Log("%c    ", LetterToChar(i));
		for (int j = 0; j < 20; ++j)
			Log("%5.1f", Mx[i][j]);
		Log("\n");
		}
	Log("\n");
	}

static unsigned MxCharToLetter(char c)
	{
	for (unsigned Letter = 0; Letter < HeadingCount; ++Letter)
		if (Heading[Letter] == c)
			return Letter;
	Quit("Letter '%c' has no heading", c);
	return 0;
	}

PTR_SCOREMATRIX ReadMx(TextFile &File)
	{
// Find column headers
	char Line[MAX_LINE];
	for (;;)
		{
		bool EndOfFile = File.GetLine(Line, sizeof(Line));
		if (EndOfFile)
			Quit("Premature EOF in matrix file");

		if (Line[0] == '#')
			continue;
		else if (Line[0] == ' ')
			break;
		else
			Quit("Invalid line in matrix file: '%s'", Line);
		}

// Read column headers
	HeadingCount = 0;
	for (char *p = Line; *p; ++p)
		{
		char c = *p;
		if (!isspace(c))
			Heading[HeadingCount++] = c;
		}

	if (HeadingCount > 0 && Heading[HeadingCount-1] == '*')
		--HeadingCount;

	if (HeadingCount < 20)
		Quit("Error in matrix file: < 20 headers, line='%s'", Line);

#if TRACE
	{
	Log("ReadMx\n");
	Log("%d headings: ", HeadingCount);
	for (unsigned i = 0; i < HeadingCount; ++i)
		Log("%c", Heading[i]);
	Log("\n");
	}
#endif

// Zero out matrix
	for (int i = 0; i < MAX_ALPHA; ++i)
		for (int j = 0; j < MAX_ALPHA; ++j)
			Mx[i][j] = 0.0;

// Read data lines
	for (unsigned RowIndex = 0; RowIndex < HeadingCount; ++RowIndex)
		{
		bool EndOfFile = File.GetTrimLine(Line, sizeof(Line));
		if (EndOfFile)
			Quit("Premature EOF in matrix file");
#if	TRACE
		Log("Line=%s\n", Line);
#endif
		if (Line[0] == '#')
			continue;

		char c = Line[0];
#if	TRACE
		Log("Row char=%c\n", c);
#endif
		if (!IsResidueChar(c))
			continue;
		unsigned RowLetter = CharToLetter(c);
		if (RowLetter >= 20)
			continue;
#if	TRACE
		Log("Row letter = %u\n", RowLetter);
#endif

		char *p = Line + 1;
		char *maxp = p + strlen(Line);
		for (unsigned Col = 0; Col < HeadingCount - 1; ++Col)
			{
			if (p >= maxp)
				Quit("Too few fields in line of matrix file: '%s'", Line);
			while (isspace(*p))
				++p;
			char *Value = p;
			while (!isspace(*p))
				++p;
			float v = (float) atof(Value);
			char HeaderChar = Heading[Col];
			if (IsResidueChar(HeaderChar))
				{
				unsigned ColLetter = CharToLetter(HeaderChar);
				if (ColLetter >= 20)
					continue;
				Mx[RowLetter][ColLetter] = v;
				}
			p += 1;
			}
		}

// Sanity check for symmetry
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < i; ++j)
			{
			if (Mx[i][j] != Mx[j][i])
				{
				Warning("Matrix is not symmetrical, %c->%c=%g, %c->%c=%g",
				  CharToLetter(i),
				  CharToLetter(j),
				  Mx[i][j],
				  CharToLetter(j),
				  CharToLetter(i),
				  Mx[j][i]);
				goto ExitLoop;
				}
			}
ExitLoop:;

	if (g_bVerbose)
		LogMx();

	return &Mx;
	}
