#ifndef	alpha_h
#define	alpha_h

enum ALPHA
	{
	ALPHA_Undefined,
	ALPHA_Nucleo,
	ALPHA_Amino
	};

bool StrHasAmino(const char *Str);
bool StrHasGap(const char *Str);
void ClearInvalidLetterWarning();
void InvalidLetterWarning(char c, char w);
void ReportInvalidLetters();

extern unsigned g_CharToLetter[];
extern unsigned g_CharToLetterEx[];

extern char g_LetterToChar[];
extern char g_LetterExToChar[];

extern char g_UnalignChar[];
extern char g_AlignChar[];

extern bool g_IsWildcardChar[];
extern bool g_IsResidueChar[];

#define CharToLetter(c)		(g_CharToLetter[(unsigned char) (c)])
#define CharToLetterEx(c)	(g_CharToLetterEx[(unsigned char) (c)])

#define LetterToChar(u)		(g_LetterToChar[u])
#define LetterExToChar(u)	(g_LetterExToChar[u])

#define IsResidueChar(c)	(g_IsResidueChar[(unsigned char) (c)])
#define IsGapChar(c)		('-' == (c) || '.' == (c))
#define IsWildcardChar(c)	(g_IsWildcardChar[(unsigned char) (c)])

#define AlignChar(c)		(g_AlignChar[(unsigned char) (c)])
#define UnalignChar(c)		(g_UnalignChar[(unsigned char) (c)])

// AX=Amino alphabet with eXtensions (B, Z and X)
enum AX
	{
	AX_A,
	AX_C,
	AX_D,
	AX_E,
	AX_F,
	AX_G,
	AX_H,
	AX_I,
	AX_K,
	AX_L,
	AX_M,
	AX_N,
	AX_P,
	AX_Q,
	AX_R,
	AX_S,
	AX_T,
	AX_V,
	AX_W,
	AX_Y,

	AX_X,	// Any

	AX_B,	// D or N
	AX_Z,	// E or Q

	AX_GAP,
	};
const unsigned AX_COUNT = AX_GAP + 1;

// NX=Nucleotide alphabet with extensions
enum NX
	{
	NX_A,
	NX_C,
	NX_G,
	NX_T,
	NX_U = NX_T,

    NX_M, // AC
    NX_R, // AG
    NX_W, // AT
    NX_S, // CG
    NX_Y, // CT
    NX_K, // GT
    NX_V, // ACG
    NX_H, // ACT
    NX_D, // AGT
    NX_B, // CGT
    NX_X, // GATC
    NX_N, // GATC
	NX_GAP
	};
const unsigned NX_COUNT = NX_GAP + 1;

const unsigned MAX_ALPHA = 20;
const unsigned MAX_ALPHA_EX = AX_COUNT;
const unsigned MAX_CHAR = 256;

extern ALPHA g_Alpha;
extern unsigned g_AlphaSize;

void SetAlpha(ALPHA Alpha);
char GetWildcardChar();
bool IsNucleo(char c);
bool IsDNA(char c);
bool IsRNA(char c);

static inline bool isgap(char c)
	{
	return c == '-' || c == '.';
	}

extern byte g_CharToLetterNucleo[256];
extern byte g_CharToLetterAmino[256];

#endif	// alpha_h
