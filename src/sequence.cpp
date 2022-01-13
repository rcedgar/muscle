#include "myutils.h"
#include "sequence.h"

#if SEQ_TRACE
#include "locallock.h"
static vector<const char *> g_Files;
static vector<int> g_Lines;
static vector<bool> g_Deleted;
static vector<const Sequence *> g_SeqPtrs;

void Sequence::AllocReport(const string &Msg)
	{
	Log("\nSequence::AllocReport(%s)\n", Msg.c_str());
	uint64 TotalBytes = 0;
	const uint N = SIZE(g_Files);
	asserta(SIZE(g_Lines) == N);
	asserta(SIZE(g_Deleted) == N);
	asserta(SIZE(g_SeqPtrs) == N);
	uint DeletedCount = 0;
	uint NotDeletedCount = 0;
	Log("\n");
	for (uint i = 0; i < N; ++i)
		{
		if (g_Deleted[i])
			continue;
		++NotDeletedCount;
		double Bytes = double(g_SeqPtrs.size());
		TotalBytes += g_SeqPtrs.capacity();
		Log("[%7u]  %10.0f bytes  %s:%u\n",
		  i, Bytes, g_Files[i], g_Lines[i]);
		}
	Log("Seqs %u / %u not freed, bytes %s\n", 
	  NotDeletedCount, N, MemBytesToStr(TotalBytes));
	}

void Sequence::AssertId() const
	{
	asserta(m_Id < SIZE(g_Files));
	}

Sequence *Sequence::_NewSequence(const char *File, int Line)
	{
	Sequence *Seq = new Sequence;
	Lock();
	uint Id = SIZE(g_Files);
	Seq->m_Id = Id;
	g_Files.push_back(File);
	g_Lines.push_back(Line);
	g_Deleted.push_back(false);
	g_SeqPtrs.push_back(Seq);
	Unlock();
	return Seq;
	}

void Sequence::_DeleteSequence(const Sequence *s,
  const char *File, int Line)
	{
	uint Id = s->m_Id;
	asserta(Id < SIZE(g_Deleted));
	asserta(!g_Deleted[Id]);
	delete s;
	g_Deleted[Id] = true;
	}

#else

Sequence *Sequence::_NewSequence()
	{
	Sequence *Seq = new Sequence;
	return Seq;
	}

void Sequence::_DeleteSequence(const Sequence *Seq)
	{
	delete Seq;
	}

#endif

void Sequence::Create(const vector<char> *a_data, string a_label, uint GSI, uint SMI)
	{
	m_CharVec = *a_data;
	m_Label = a_label;
	m_GSI = GSI;
	m_SMI = SMI;
	}

bool Sequence::FromFileBuffer(FileBuffer& infile, bool stripGaps)
	{
	if (infile.eof())
		return false;

	m_GSI = UINT_MAX;
	m_SMI = UINT_MAX;
	m_Label = "~";

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
	m_CharVec.push_back('@');

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

		if (stripGaps && ch == '-')
			continue;

		m_CharVec.push_back(ch);
		}
	return true;
	}

void Sequence::WriteMFA(FILE *f) const
	{
	const vector<char> &v = m_CharVec;
	const int L = GetLength();
	byte *Seq = myalloc(byte, L);
	for (int i = 0; i < L; ++i)
		Seq[i] = v[i+1];
	SeqToFasta(f, Seq, (uint) L, m_Label.c_str());
	}

Sequence* Sequence::Clone() const
	{
	Sequence* ret = NewSequence();
	asserta(ret);

	ret->m_Label = m_Label;
	ret->m_CharVec = m_CharVec;
	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;

	return ret;
	}

/////////////////////////////////////////////////////////////////
// Sequence::AddGaps()
//
// Given an vector<char> containing the skeleton for an
// alignment and the identity of the current character, this
// routine will create a new sequence with all necesssary gaps added.
// For instance,
//    alignment = "XXXBBYYYBBYYXX"
//    id = 'X'
// will perform the transformation
//    "ATGCAGTCA" --> "ATGCC---GT--CA"
//                    (XXXBBYYYBBYYXX)
/////////////////////////////////////////////////////////////////
Sequence* Sequence::AddGaps(const vector<char>* alignment, char id) const
	{
	Sequence* ret = NewSequence();
	assert(ret);

	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;

	ret->m_Label = m_Label;
	ret->m_CharVec.clear();
	ret->m_CharVec.push_back('@');

	vector<char>::const_iterator dataIter = m_CharVec.begin() + 1;
	for (vector<char>::const_iterator iter = alignment->begin();
	  iter != alignment->end(); ++iter)
		{
		if (*iter == 'B' || *iter == id)
			{
			ret->m_CharVec.push_back(*dataIter);
			++dataIter;
			}
		else
			ret->m_CharVec.push_back('-');
		}

	return ret;
	}

Sequence* Sequence::AddGapsPath(const string &Path, char id) const
	{
	Sequence* ret = NewSequence();
	assert(ret);

	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;

	ret->m_Label = m_Label;
	ret->m_CharVec.clear();
	ret->m_CharVec.push_back('@');

	vector<char>::const_iterator dataIter = m_CharVec.begin() + 1;
	for (string::const_iterator iter = Path.begin();
	  iter != Path.end(); ++iter)
		{
		if (*iter == 'B' || *iter == id)
			{
			ret->m_CharVec.push_back(*dataIter);
			++dataIter;
			}
		else
			ret->m_CharVec.push_back('-');
		}

	return ret;
	}

// Returns vector containing 1-based col indexes, e.g. if m_CharVec
// is "ATGCC---GT--CA" vector is set to {1,2,3,4,5,9,10,13,14}.
void Sequence::GetPosToCol_OneBased(vector<uint> &PosToCol) const
	{
	PosToCol.clear();
	PosToCol.push_back(UINT_MAX);
	uint L = GetLength();
	for (uint i = 1; i <= L; i++)
		{
		if (m_CharVec[i] != '-')
			PosToCol.push_back(i);
		}
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

Sequence *Sequence::DeleteGaps() const
	{
	Sequence* ret = NewSequence();
	asserta(ret);

	ret->m_Label = m_Label;
	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;
	ret->m_CharVec.clear();
	asserta(!m_CharVec.empty() && m_CharVec[0] == '@');
	int L = GetLength();
	ret->m_CharVec.push_back('@');
	for (int i = 1; i <= L; ++i)
		{
		char c = m_CharVec[i];
		if (c != '-')
			ret->m_CharVec.push_back(c);
		}
	return ret;
	}

//void Sequence::Copy(const Sequence &rhs)
//	{
//	m_Label = rhs.m_Label;
//	m_SMI = rhs.m_SMI;
//	m_GSI = rhs.m_GSI;
//	m_CharVec = rhs.m_CharVec;
//	}

void Sequence::FromString(const string &Label, const string &Seq)
	{
	m_GSI = UINT_MAX;
	m_SMI = UINT_MAX;

	m_Label = Label;
	m_CharVec.clear();
	int L = (int) SIZE(Seq);
	m_CharVec.push_back('@');
	for (int i = 0; i < L; ++i)
		{
		char c = Seq[i];
		m_CharVec.push_back(c);
		}
	}

void Sequence::LogMe() const
	{
	Log("\n");
	uint L = GetLength();
	Log("Sequence(%p) length %u  >%s\n",
	  this, L, m_Label.c_str());
	for (uint i = 1; i <= L; ++i)
		Log("%c", m_CharVec[i]);
	Log("\n");
	}
