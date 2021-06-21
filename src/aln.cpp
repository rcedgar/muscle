#include "muscle.h"
#include <stdio.h>
#include <ctype.h>
#include "msa.h"
#include "textfile.h"

const unsigned uCharsPerLine = 60;
const int MIN_NAME = 10;
const int MAX_NAME = 32;

static char GetAlnConsensusChar(const MSA &a, unsigned uColIndex);

void MSA::ToAlnFile(TextFile &File) const
	{
	if (g_bClwStrict)
		File.PutString("CLUSTAL W (1.81) multiple sequence alignment\n");
	else
		{
		File.PutString("MUSCLE ("
		  SHORT_VERSION ")"
		  " multiple sequence alignment\n");
		File.PutString("\n");
		}

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

			File.PutFormat("%s      ", Name);
			for (unsigned uColIndex = uStartColIndex; uColIndex <= uEndColIndex;
			  ++uColIndex)
				{
				const char c = GetChar(uSeqIndex, uColIndex);
				File.PutFormat("%c", toupper(c));
				}
			File.PutString("\n");
			}

		memset(Name, ' ', MAX_NAME);
		Name[iLongestNameLength] = 0;
		File.PutFormat("%s      ", Name);
		for (unsigned uColIndex = uStartColIndex; uColIndex <= uEndColIndex;
		  ++uColIndex)
			{
			const char c = GetAlnConsensusChar(*this, uColIndex);
			File.PutChar(c);
			}
		File.PutString("\n");
		}
	}

static char GetAlnConsensusChar(const MSA &a, unsigned uColIndex)
	{
	const unsigned uSeqCount = a.GetSeqCount();
	unsigned BitMap = 0;
	unsigned Count = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uLetter = a.GetLetterEx(uSeqIndex, uColIndex);
		assert(uLetter < 32);
		unsigned Bit = (1 << uLetter);
		if (!(BitMap & Bit))
			++Count;
		BitMap |= Bit;
		}

//	'*' indicates positions which have a single, fully conserved residue
	if (1 == Count)
		return '*';

	if (ALPHA_Amino != g_Alpha)
		return ' ';

#define B(a)	(1 << AX_##a)
#define S2(a, b)		S(B(a) | B(b))
#define S3(a, b, c)		S(B(a) | B(b) | B(c))
#define S4(a, b, c, d)	S(B(a) | B(b) | B(c) | B(d))
#define S(w)	if (0 == (BitMap & ~(w)) && (BitMap & (w)) != 0) return ':';

#define W3(a, b, c)				W(B(a) | B(b) | B(c))
#define W4(a, b, c, d)			W(B(a) | B(b) | B(c) | B(d))
#define W5(a, b, c, d, e)		W(B(a) | B(b) | B(c) | B(d) | B(e))
#define W6(a, b, c, d, e, f)	W(B(a) | B(b) | B(c) | B(d) | B(e) | B(f))
#define W(w)	if (0 == (BitMap & ~(w)) && (BitMap & (w)) != 0) return '.';

//	':' indicates that one of the following 'strong'
// groups is fully conserved
//                 STA
//                 NEQK
//                 NHQK
//                 NDEQ
//                 QHRK
//                 MILV
//                 MILF
//                 HY
//                 FYW
//
	S3(S, T, A)
	S4(N, E, Q, K)
	S4(N, H, Q, K)
	S4(N, D, E, Q)
	S4(M, I, L, V)
	S4(M, I, L, F)
	S2(H, Y)
	S3(F, Y, W)

//	'.' indicates that one of the following 'weaker' 
// groups is fully conserved
//                 CSA
//                 ATV
//                 SAG
//                 STNK
//                 STPA
//                 SGND
//                 SNDEQK
//                 NDEQHK
//                 NEQHRK
//                 FVLIM
//                 HFY
	W3(C, S, A)
	W3(A, T, V)
	W3(S, A, G)
	W4(S, T, N, K)
	W4(S, T, P, A)
	W4(S, G, N, D)
	W6(S, N, D, E, Q, K)
	W6(N, W, Q, H, R, K)
	W5(F, V, L, I, M)
	W3(H, F, Y)

	return ' ';
	}
