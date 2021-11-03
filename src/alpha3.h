#ifndef alpha3_h
#define alpha3_h

#include <limits.h>
#include <string>

using namespace std;

const byte INVALID_LETTER = 0xff;
const byte INVALID_CHAR = '?';
const unsigned BAD_WORD = UINT_MAX;

extern byte g_AminoAcidChars[];
extern byte g_CharToLetterAmino[];
extern byte g_CharToLetterAminoStop[];
extern byte g_CharToLetterAminoGap[];
extern byte g_LetterToCharAmino[];
extern byte g_LetterToCharAminoGap[];
extern byte g_CharToLetterNucleo[];
extern byte g_CharToLetterNucleoGap[];
extern byte g_CharToLetterNucleoMasked[];
extern byte g_LetterToCharNucleo[];
extern byte g_LetterToCharNucleoGap[];
extern byte g_CodonWordToAminoLetter[];
extern byte g_CodonWordToAminoChar[];
extern byte g_CharToCompChar[];
extern byte g_CharToCompLetter[];
extern byte g_IUPAC_PairCharToChar1[256];
extern byte g_IUPAC_PairCharToChar2[256];
extern byte g_IUPAC_PairCharToCharCase[256];
extern byte g_CharToLetterSEB8[256];
extern bool **g_MatchMxNucleo;
extern bool **g_MatchMxAmino;

extern bool g_IsAminoChar[];
extern bool g_IsNucleoChar[];
extern bool g_IsACGTU[];
extern bool g_IsSeqChar[];

extern float g_AminoFreqs[];

extern unsigned g_CharToLetterRed[];
extern byte g_LetterToCharRed[];
extern unsigned g_RedAlphaSize;

void LogRedAlphaRed();
void ReadRedAlphaFromFile(const string &FileName);
byte GetAminoCharFrom3NucChars(byte c1, byte c2,
  byte c3);

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo);
const char *WordToStrNucleo(unsigned Word, unsigned WordLength);
const char *WordToStrAmino(unsigned Word, unsigned WordLength);
const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str);

static inline bool isgap(byte c)
	{
	return c == '-' || c == '.';
	}

void InitAlpha();
byte IUPAC_Pair(byte CharOrWildcard1, byte CharOrWildcard2);

#endif // alpha3_h
