#ifndef alpha3_h
#define alpha3_h

#include <limits.h>
#include <string>

using namespace std;

const unsigned INVALID_LETTER = 0xff;
const unsigned char INVALID_CHAR = '?';
const unsigned BAD_WORD = UINT_MAX;

extern unsigned char g_AminoAcidChars[];
extern unsigned char g_CharToLetterAmino[];
extern unsigned char g_CharToLetterAminoStop[];
extern unsigned char g_CharToLetterAminoGap[];
extern unsigned char g_LetterToCharAmino[];
extern unsigned char g_LetterToCharAminoGap[];
extern unsigned char g_CharToLetterNucleo[];
extern unsigned char g_CharToLetterNucleoGap[];
extern unsigned char g_CharToLetterNucleoMasked[];
extern unsigned char g_LetterToCharNucleo[];
extern unsigned char g_LetterToCharNucleoGap[];
extern unsigned char g_CodonWordToAminoLetter[];
extern unsigned char g_CodonWordToAminoChar[];
extern unsigned char g_CharToCompChar[];
extern unsigned char g_CharToCompLetter[];
extern unsigned char g_IUPAC_PairCharToChar1[256];
extern unsigned char g_IUPAC_PairCharToChar2[256];
extern unsigned char g_IUPAC_PairCharToCharCase[256];
extern bool **g_MatchMxNucleo;
extern bool **g_MatchMxAmino;

extern bool g_IsAminoChar[];
extern bool g_IsNucleoChar[];
extern bool g_IsACGTU[];
extern bool g_IsSeqChar[];

extern float g_AminoFreqs[];

extern unsigned g_CharToLetterRed[];
extern unsigned char g_LetterToCharRed[];
extern unsigned g_RedAlphaSize;

void LogRedAlphaRed();
void ReadRedAlphaFromFile(const string &FileName);
unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3);

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo);
const char *WordToStrNucleo(unsigned Word, unsigned WordLength);
const char *WordToStrAmino(unsigned Word, unsigned WordLength);
const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str);

static inline bool isgap(unsigned char c)
	{
	return c == '-' || c == '.';
	}

void InitAlpha();
unsigned char IUPAC_Pair(unsigned char CharOrWildcard1, unsigned char CharOrWildcard2);

#endif // alpha3_h
