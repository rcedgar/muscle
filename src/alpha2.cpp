#include "muscle.h"

//void SetSubstMx(bool IsNucleo);

ALPHA g_Alpha = ALPHA_Undefined;
unsigned g_AlphaSize = 0;

byte g_CharToLetter[256];
byte g_LetterToChar[256];

void SetAlphaLC(bool IsNucleo)
	{
	if (IsNucleo)
		{
		if (g_Alpha == ALPHA_Nucleo)
			return;
		asserta(g_Alpha == ALPHA_Undefined);
		g_Alpha = ALPHA_Nucleo;
		memcpy(g_CharToLetter, g_CharToLetterNucleo, 256);
		memcpy(g_LetterToChar, g_LetterToCharNucleo, 256);
		g_AlphaSize = 4;
		}
	else
		{
		if (g_Alpha == ALPHA_Amino)
			return;
		asserta(g_Alpha == ALPHA_Undefined);
		g_Alpha = ALPHA_Amino;
		memcpy(g_CharToLetter, g_CharToLetterAmino, 256);
		memcpy(g_LetterToChar, g_LetterToCharAmino, 256);
		g_AlphaSize = 20;
		}
	}

void SetAlpha(ALPHA Alpha)
	{
	switch (Alpha)
		{
	case ALPHA_Amino:
		SetAlphaLC(false);
		break;

	case ALPHA_Nucleo:
		SetAlphaLC(true);
		break;

	default:
		Die("Invalid Alpha=%d", Alpha);
		}
	}

void SetAlphab(bool IsNucleo)
	{
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);
	}
