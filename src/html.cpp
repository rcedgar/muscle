#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

const unsigned uCharsPerLine = 60;
const int MIN_NAME = 10;
const int MAX_NAME = 32;

extern void AssignColors(const MSA &a, int **Colors);

static int **MakeColors(const MSA &a)
	{
	const unsigned uSeqCount = a.GetSeqCount();
	const unsigned uColCount = a.GetColCount();

	int **Colors = new int *[uSeqCount];
	for (unsigned i = 0; i < uSeqCount; ++i)
		{
		Colors[i] = new int[uColCount];
		memset(Colors[i], 0, uColCount*sizeof(int));
		}
	AssignColors(a, Colors);
	return Colors;
	}

static void ChangeColor(TextFile &File, int From, int To)
	{
	if (From == To)
		return;

#define	COLOR_WHITE		"FFFFFF"
#define	COLOR_GRAY		"C0C0C0"
#define	COLOR_BLACK		"000000"
#define COLOR_RED		"FF0000"
#define COLOR_GREEN		"00FF00"
#define COLOR_BLUE		"5590FF"
#define COLOR_LIGHTBLUE	"77FFFF"

#define X(c)	File.PutString("</SPAN><SPAN STYLE=\"background-color:#" c "\">");
	switch (To)
		{
	case 0:
		X(COLOR_WHITE)
		break;
	case 1:
		X(COLOR_GRAY)
		break;
	case 2:
		X(COLOR_BLUE)
		break;
	case 3:
		X(COLOR_LIGHTBLUE)
		break;
		}
	}

#define COLOR_WINDOW "FFEEE0"

void MSA::ToHTMLFile(TextFile &File) const
	{
	File.PutString("<HTML>\n");
	File.PutString("<BODY BGCOLOR=\"#" COLOR_WINDOW "\">\n");
	File.PutString("<PRE>");

	int **Colors = MakeColors(*this);

	int iLongestNameLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char *ptrName = GetSeqName(uSeqIndex);
		const char *ptrBlank = strchr(ptrName, ' ');
		int iLength;
		if (0 != ptrBlank)
			iLength = (int) (ptrBlank - ptrName);
		else
			iLength = (int) strlen(ptrName);
		if (iLength > iLongestNameLength)
			iLongestNameLength = iLength;
		}
	if (iLongestNameLength > MAX_NAME)
		iLongestNameLength = MAX_NAME;
	if (iLongestNameLength < MIN_NAME)
		iLongestNameLength = MIN_NAME;

	unsigned uLineCount = (GetColCount() - 1)/uCharsPerLine + 1;
	int CurrentColor = -1;
	for (unsigned uLineIndex = 0; uLineIndex < uLineCount; ++uLineIndex)
		{
		File.PutString("\n");
		unsigned uStartColIndex = uLineIndex*uCharsPerLine;
		unsigned uEndColIndex = uStartColIndex + uCharsPerLine - 1;
		if (uEndColIndex >= GetColCount())
			uEndColIndex = GetColCount() - 1;
		char Name[MAX_NAME+1];
		for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
			{
			const char *ptrName = GetSeqName(uSeqIndex);
			const char *ptrBlank = strchr(ptrName, ' ');
			int iLength;
			if (0 != ptrBlank)
				iLength = (int) (ptrBlank - ptrName);
			else
				iLength = (int) strlen(ptrName);
			if (iLength > MAX_NAME)
				iLength = MAX_NAME;
			memset(Name, ' ', MAX_NAME);
			memcpy(Name, ptrName, iLength);
			Name[iLongestNameLength] = 0;

//			File.PutString("<FONT COLOR=\"#000000\">");
			CurrentColor = -1;
			File.PutString("<SPAN STYLE=\"background-color:#" COLOR_WINDOW "\">");
			File.PutFormat("%s      ", Name);
			File.PutString("<SPAN STYLE=\"background-color:#FFFFFF\">");
			for (unsigned uColIndex = uStartColIndex; uColIndex <= uEndColIndex;
			  ++uColIndex)
				{
				const int Color = Colors[uSeqIndex][uColIndex];
				ChangeColor(File, CurrentColor, Color);
				CurrentColor = Color;
				const char c = GetChar(uSeqIndex, uColIndex);
				if (Color == 0)
					File.PutFormat("%c", tolower(c));
				else
					File.PutFormat("%c", toupper(c));
				}
			File.PutString("\n");
			}
		}
	File.PutString("</SPAN>\n");
	File.PutString("</PRE>\n");
	File.PutString("</BODY>\n");
	File.PutString("</HTML>\n");
	}
