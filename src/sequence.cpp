#include "myutils.h"
#include "sequence.h"
#include "alpha.h"

static uint g_NewCount;
static uint g_DeleteCount;

void Sequence::LogNewDeleteCounts()
	{
	Log("Sequence::LogNewDeleteCounts new=%u, delete=%u\n",
	  g_NewCount, g_DeleteCount);
	}

Sequence *Sequence::_NewSequence()
	{
	Sequence *Seq = new Sequence;
#pragma omp critical
	++g_NewCount;

	return Seq;
	}

void Sequence::_DeleteSequence(const Sequence *Seq)
	{
#pragma omp critical
	++g_DeleteCount;

	delete Seq;
	}

//void Sequence::Create(const vector<char> *a_data, string a_label)
//	{
//	m_CharVec = *a_data;
//	m_Label = a_label;
//	}

bool Sequence::FromFileBuffer(FileBuffer& infile, bool stripGaps)
	{
	if (infile.eof())
		return false;

	m_Label.clear();

// Skip blank lines
	for (;;)
		{
		if (infile.eof())
			{
			if (m_Label.empty())
				return false;
			asserta(false);
			}
		infile.GetLine(m_Label);
		if (m_Label.length() > 0)
			break;
		}

	if (m_Label[0] != '>')
		Die("Expected '>' in FASTA, got '%s'", m_Label.c_str());

// Remove leading ">"
	m_Label = m_Label.substr(1);

	if (opt(accs))
		{
		string Acc;
		GetAccFromLabel(m_Label, Acc);
		m_Label = Acc;
		}

	m_CharVec.clear();

	char ch;
	while (infile.Get(ch))
		{
		if (ch == '>')
			{
			infile.UnGet();
			break;
			}

		if (isspace(ch)) 
			continue;

		if (stripGaps && (ch == '-' || ch == '.'))
			continue;
		if (stripGaps)
			ch = toupper(ch);

		m_CharVec.push_back(ch);
		}
	return true;
	}

void Sequence::WriteMFA(FILE *f) const
	{
	const byte *Seq = GetBytePtr();
	const int L = GetLength();
	SeqToFasta(f, Seq, (uint) L, m_Label.c_str());
	}

Sequence* Sequence::Clone() const
	{
	Sequence* ret = NewSequence();
	asserta(ret);

	ret->m_Label = m_Label;
	ret->m_CharVec = m_CharVec;
	//ret->m_GSI = m_GSI;
	//ret->m_SMI = m_SMI;

	return ret;
	}

Sequence* Sequence::AddGapsPath(const string &Path, char id) const
	{
	Sequence* ret = NewSequence();
	assert(ret);

	//ret->m_GSI = m_GSI;
	//ret->m_SMI = m_SMI;

	ret->m_Label = m_Label;
	ret->m_CharVec.clear();

	vector<char>::const_iterator dataIter = m_CharVec.begin();
	for (string::const_iterator iter = Path.begin();
	  iter != Path.end(); ++iter)
		{
		if (*iter == 'M' || *iter == 'B' || *iter == id)
			{
			ret->m_CharVec.push_back(*dataIter);
			++dataIter;
			}
		else
			ret->m_CharVec.push_back('-');
		}

	return ret;
	}

// Returns vector containing 0-based col indexes, e.g.
// "ATGCC---GT--CA" PosToCol={0,1,2,3,4,8,9,12,13}.
void Sequence::GetPosToCol(vector<uint> &PosToCol) const
	{
	PosToCol.clear();
	const uint ColCount = GetLength();
	const byte *ByteSeq = GetBytePtr();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (ByteSeq[Col] != '-')
			PosToCol.push_back(Col);
		}
	}

void Sequence::GetSeqAsString(string &Seq) const
	{
	Seq.clear();
	const uint L = GetLength();
	const char *CharSeq = GetCharPtr();
	for (uint i = 0; i < L; ++i)
		Seq += CharSeq[i];
	}

void Sequence::GetColToPos(vector<uint> &ColToPos) const
	{
	ColToPos.clear();
	const uint ColCount = GetLength();
	const byte *ByteSeq = GetBytePtr();
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (ByteSeq[Col] == '-')
			ColToPos.push_back(UINT_MAX);
		else
			ColToPos.push_back(Pos++);
		}
	}

Sequence *Sequence::CopyDeleteGaps() const
	{
	Sequence* ret = NewSequence();
	asserta(ret);

	ret->m_Label = m_Label;
	ret->m_CharVec.clear();
	int L = GetLength();
	for (int i = 0; i < L; ++i)
		{
		char c = m_CharVec[i];
		if (c != '-')
			ret->m_CharVec.push_back(c);
		}
	return ret;
	}

void Sequence::FromString(const string &Label, const string &Seq)
	{
	m_Label = Label;
	m_CharVec.clear();
	int L = (int) SIZE(Seq);
	for (int i = 0; i < L; ++i)
		{
		char c = Seq[i];
		m_CharVec.push_back(c);
		}
	}

void Sequence::LogMe() const
	{
	uint L = GetLength();
	for (uint i = 0; i < L; ++i)
		Log("%c", m_CharVec[i]);
	Log("  >%s (%u)\n", m_Label.c_str(), L);
	}
