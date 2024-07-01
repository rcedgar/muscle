#ifndef	alpha_h
#define	alpha_h

enum ALPHA
	{
	ALPHA_Undefined,
	ALPHA_Nucleo,
	ALPHA_Amino
	};

void ClearInvalidLetterWarning();
void InvalidLetterWarning(char c, char w);
void ReportInvalidLetters();

extern byte g_CharToLetterEx[];

extern byte g_LetterToCharAmino[];
extern byte g_LetterToChar[];
extern byte g_LetterExToChar[];

#define CharToLetter(c)		(g_CharToLetter[(unsigned char) (c)])
#define CharToLetterEx(c)	(g_CharToLetterEx[(unsigned char) (c)])

#define LetterToChar(u)		(g_LetterToChar[u])

#define IsResidueChar(c)	(g_IsResidueChar[(unsigned char) (c)])
#define IsGapChar(c)		('-' == (c) || '.' == (c))
#define IsWildcardChar(c)	(g_IsWildcardChar[(unsigned char) (c)])

#define AlignChar(c)		(g_AlignChar[(unsigned char) (c)])
#define UnalignChar(c)		(g_UnalignChar[(unsigned char) (c)])

const unsigned MAX_CHAR = 256;

extern ALPHA g_Alpha;
extern unsigned g_AlphaSize;

void SetAlphaLC(bool IsNucleo);
void SetAlpha(ALPHA Alpha);
void SetAlphab(bool IsNucleo);
char GetWildcardChar();
bool CharIsNucleo(char c);
bool CharIsDNA(char c);
bool CharIsRNA(char c);

static inline bool isgap(char c)
	{
	return c == '-' || c == '.';
	}

extern byte g_CharToLetterNucleo[256];
extern byte g_CharToLetterAmino[256];
extern byte g_CharToLetter[256];

#endif	// alpha_h
