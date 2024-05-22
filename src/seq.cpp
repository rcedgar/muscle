#include "muscle.h"
#include "seq.h"
#include "textfile.h"
#include "msa.h"
//#include <ctype.h>

const size_t MAX_FASTA_LINE = 16000;

void Seq::SetName(const char *ptrName)
	{
	delete[] m_ptrName;
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	}

void Seq::ToFASTAFile(TextFile &File) const
	{
	File.PutFormat(">%s\n", m_ptrName);
	unsigned uColCount = Length();
	for (unsigned n = 0; n < uColCount; ++n)
		{
		if (n > 0 && n%60 == 0)
			File.PutString("\n");
		File.PutChar(at(n));
		}
	File.PutString("\n");
	}

// Return true on end-of-file
bool Seq::FromFASTAFile(TextFile &File)
	{
	Clear();

	char szLine[MAX_FASTA_LINE];
	bool bEof = File.GetLine(szLine, sizeof(szLine));
	if (bEof)
		return true;
	if ('>' != szLine[0])
		Die("Expecting '>' in FASTA file %s line %u",
		  File.GetFileName(), File.GetLineNr());

	size_t n = strlen(szLine);
	if (1 == n)
		Die("Missing annotation following '>' in FASTA file %s line %u",
		  File.GetFileName(), File.GetLineNr());

	m_ptrName = new char[n];
	strcpy(m_ptrName, szLine + 1);

	TEXTFILEPOS Pos = File.GetPos();
	for (;;)
		{
		bEof = File.GetLine(szLine, sizeof(szLine));
		if (bEof)
			{
			if (0 == size())
				{
				Die("Empty sequence in FASTA file %s line %u",
				  File.GetFileName(), File.GetLineNr());
				return true;
				}
			return false;
			}
		if ('>' == szLine[0])
			{
			if (0 == size())
				Die("Empty sequence in FASTA file %s line %u",
				  File.GetFileName(), File.GetLineNr());
		// Rewind to beginning of this line, it's the start of the
		// next sequence.
			File.SetPos(Pos);
			return false;
			}
		const char *ptrChar = szLine;
		while (char c = *ptrChar++)
			{
			if (isspace(c))
				continue;
			if (IsGapChar(c))
				continue;
			if (!isalpha(c))
				{
				static bool WarningDone = false;
				if (!WarningDone)
					{
					if (isprint(c))
						Warning("Invalid character'%c' in FASTA file %s line %d",
						  c, File.GetFileName(), File.GetLineNr());
					else
						Warning("Invalid byte hex %02x in FASTA file %s line %d",
						  (unsigned char) c, File.GetFileName(), File.GetLineNr());
					}
				continue;
				}
			c = toupper(c);
			push_back(c);
			}
		Pos = File.GetPos();
		}
	}

void Seq::ExtractUngapped(MSA &msa) const
	{
	msa.Clear();
	unsigned uColCount = Length();
	msa.SetSize(1, 1);
	unsigned uUngappedPos = 0;
	for (unsigned n = 0; n < uColCount; ++n)
		{
		char c = at(n);
		if (!IsGapChar(c))
			msa.SetChar(0, uUngappedPos++, c);
		}
	msa.SetSeqName(0, m_ptrName);
	}

void Seq::Copy(const Seq &rhs)
	{
	clear();
	const unsigned uLength = rhs.Length();
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(rhs.at(uColIndex));
	const char *ptrName = rhs.GetName();
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	SetId(rhs.GetId());
	}

void Seq::CopyReversed(const Seq &rhs)
	{
	clear();
	const unsigned uLength = rhs.Length();
	const unsigned uBase = rhs.Length() - 1;
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(rhs.at(uBase - uColIndex));
	const char *ptrName = rhs.GetName();
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	}

void Seq::StripGaps()
	{
	for (CharVect::iterator p = begin(); p != end(); )
		{
		char c = *p;
		if (IsGapChar(c))
			erase(p);
		else
			++p;
		}
	}

void Seq::StripGapsAndWhitespace()
	{
	for (CharVect::iterator p = begin(); p != end(); )
		{
		char c = *p;
		if (isspace(c) || IsGapChar(c))
			erase(p);
		else
			++p;
		}
	}

void Seq::ToUpper()
	{
	for (CharVect::iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (islower(c))
			*p = toupper(c);
		}
	}

unsigned Seq::GetLetter(unsigned uIndex) const
	{
	assert(uIndex < Length());
	char c = operator[](uIndex);
	return CharToLetter(c);
	}

bool Seq::EqIgnoreCase(const Seq &s) const
	{
	const unsigned n = Length();
	if (n != s.Length())
		return false;
	for (unsigned i = 0; i < n; ++i)
		{
		const char c1 = at(i);
		const char c2 = s.at(i);
		if (IsGapChar(c1))
			{
			if (!IsGapChar(c2))
				return false;
			}
		else
			{
			if (toupper(c1) != toupper(c2))
				return false;
			}
		}
	return true;
	}

bool Seq::Eq(const Seq &s) const
	{
	const unsigned n = Length();
	if (n != s.Length())
		return false;
	for (unsigned i = 0; i < n; ++i)
		{
		const char c1 = at(i);
		const char c2 = s.at(i);
		if (c1 != c2)
			return false;
		}
	return true;
	}

bool Seq::EqIgnoreCaseAndGaps(const Seq &s) const
	{
	const unsigned uThisLength = Length();
	const unsigned uOtherLength = s.Length();
	
	unsigned uThisPos = 0;
	unsigned uOtherPos = 0;

	int cThis;
	int cOther;
	for (;;)
		{
		if (uThisPos == uThisLength && uOtherPos == uOtherLength)
			break;

	// Set cThis to next non-gap character in this string
	// or -1 if end-of-string.
		for (;;)
			{
			if (uThisPos == uThisLength)
				{
				cThis = -1;
				break;
				}
			else
				{
				cThis = at(uThisPos);
				++uThisPos;
				if (!IsGapChar(cThis))
					{
					cThis = toupper(cThis);
					break;
					}
				}
			}

	// Set cOther to next non-gap character in s
	// or -1 if end-of-string.
		for (;;)
			{
			if (uOtherPos == uOtherLength)
				{
				cOther = -1;
				break;
				}
			else
				{
				cOther = s.at(uOtherPos);
				++uOtherPos;
				if (!IsGapChar(cOther))
					{
					cOther = toupper(cOther);
					break;
					}
				}
			}

	// Compare characters are corresponding ungapped position
		if (cThis != cOther)
			return false;
		}
	return true;
	}

unsigned Seq::GetUngappedLength() const
	{
	unsigned uUngappedLength = 0;
	for (CharVect::const_iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (!IsGapChar(c))
			++uUngappedLength;
		}
	return uUngappedLength;
	}

void Seq::LogMe() const
	{
	Log(">%s\n", m_ptrName);
	const unsigned n = Length();
	for (unsigned i = 0; i < n; ++i)
		Log("%c", at(i));
	Log("\n");
	}

void Seq::FromString(const char *pstrSeq, const char *pstrName)
	{
	clear();
	const unsigned uLength = (unsigned) strlen(pstrSeq);
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(pstrSeq[uColIndex]);
	size_t n = strlen(pstrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, pstrName);
	}

bool Seq::HasGap() const
	{
	for (CharVect::const_iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (IsGapChar(c))
			return true;
		}
	return false;
	}

//void Seq::FixAlpha()
//	{
//	for (CharVect::iterator p = begin(); p != end(); ++p)
//		{
//		char c = *p;
//		if (!IsResidueChar(c))
//			{
//			char w = GetWildcardChar();
//			// Warning("Invalid residue '%c', replaced by '%c'", c, w);
//			InvalidLetterWarning(c, w);
//			*p = w;
//			}
//		}
//	}
