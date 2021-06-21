#include "myutils.h"
#include "sequence.h"

void Sequence::Create(vector<char>* a_data, string a_label, uint GSI, uint SMI)
	{
	assert(a_data);
	assert((*a_data)[0] == '@');

	data = a_data;
	label = a_label;
	m_GSI = GSI;
	m_SMI = SMI;
	}

bool Sequence::FromFileBuffer(FileBuffer& infile, bool stripGaps)
	{
	if (infile.eof())
		return false;

	label = "~";
	asserta(data == 0);
	asserta(m_GSI == UINT_MAX);
	asserta(m_SMI == UINT_MAX);

// Skip blank lines
	for (;;)
		{
		if (infile.eof())
			{
			if (label.empty())
				return false;
			asserta(false);
			}
		infile.GetLine(label);
		if (label.length() > 0)
			break;
		}

	if (label[0] != '>')
		Die("Expected '>' in FASTA, got '%s'", label.c_str());

// Remove leading ">"
	label = label.substr(1);

	if (opt(accs))
		{
		string Acc;
		GetAccFromLabel(label, Acc);
		label = Acc;
		}

	data = new vector<char>; assert(data);
	data->push_back('@');

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

		data->push_back(ch);
		}
	return true;
	}

// Writes the sequence to outfile in MFA format.  Uses numColumns
// columns per line.  If useIndex is set to false, then the
// label is printed as normal, but if useIndex is true, then
// ">S###" is printed where ### represents the sequence label.
void Sequence::WriteMFA(FILE *f) const
	{
	const vector<char> &v = (*data);
	const int L = GetLength();
	byte *Seq = myalloc(byte, L);
	for (int i = 0; i < L; ++i)
		Seq[i] = v[i+1];
	SeqToFasta(f, Seq, (uint) L, label.c_str());
	}

Sequence* Sequence::Clone() const
	{
	Sequence* ret = new Sequence();
	asserta(ret);

	ret->label = label;
	ret->data = new vector<char>;
	asserta(ret->data);
	*(ret->data) = *data;
	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;

	return ret;
	}

//// Extract subset of positions using ZERO-based column indexes
//// (note conflict between RCE convention of 0-based and CBD
//// convention of 1-based).
//Sequence* Sequence::MakeSubset(const vector<uint> &Cols) const
//	{
//	Sequence* ret = new Sequence();
//	asserta(ret);
//
//	const uint ColCount = SIZE(Cols);
//	asserta(ColCount > 0);
//
//	ret->label = label;
//	ret->data = new vector<char>;
//	asserta(ret->data);
//	//ret->sequenceIndex = sequenceIndex;
//	ret->inputIndex = inputIndex;
//	ret->data->push_back('@');
//
//	const char *ZeroBasedSeqData = data->data() + 1;
//	asserta(ColCount <= (uint) GetLength());
//	for (uint i = 0; i < ColCount; ++i)
//		{
//		uint Pos = Cols[i];
//		char c = ZeroBasedSeqData[Pos];
//		ret->data->push_back(c);
//		}
//
//	return ret;
//	}

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
	Sequence* ret = new Sequence();
	assert(ret);

	ret->SetGSI(m_GSI);
	ret->SetSMI(m_SMI);

	ret->label = label;
	ret->data = new vector<char>; 
	assert(ret->data != 0);
	ret->data->push_back('@');

	vector<char>::iterator dataIter = data->begin() + 1;
	for (vector<char>::const_iterator iter = alignment->begin();
	  iter != alignment->end(); ++iter)
		{
		if (*iter == 'B' || *iter == id)
			{
			ret->data->push_back(*dataIter);
			++dataIter;
			}
		else
			ret->data->push_back('-');
		}

	return ret;
	}

// Returns vector containing 1-based col indexes, e.g. if data
// is "ATGCC---GT--CA" vector is set to {1,2,3,4,5,9,10,13,14}.
void Sequence::GetPosToCol_OneBased(vector<uint> &PosToCol) const
	{
	PosToCol.clear();
	PosToCol.push_back(UINT_MAX);
	uint L = GetLength();
	for (uint i = 1; i <= L; i++)
		{
		if ((*data)[i] != '-')
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
	Sequence* ret = new Sequence();
	asserta(ret);

	ret->label = label;
	ret->m_GSI = m_GSI;
	ret->m_SMI = m_SMI;
	ret->data = new vector<char>;
	asserta(ret->data);
	asserta(!(*data).empty() && (*data)[0] == '@');
	int L = GetLength();
	(*(ret->data)).push_back('@');
	for (int i = 1; i <= L; ++i)
		{
		char c = (*data)[i];
		if (c != '-')
			(*(ret->data)).push_back(c);
		}
	return ret;
	}

void Sequence::Copy(const Sequence &rhs)
	{
	label = rhs.label;
	m_GSI = rhs.m_GSI;
	m_SMI = rhs.m_SMI;

	if (rhs.data == 0)
		data = 0;
	else
		{
		data = new vector<char>;
		asserta(data);
		const vector<char> &rhsdata = *rhs.data;
		int L = rhs.GetLength();
		for (int i = 0; i <= L; ++i)
			{
			char c = (*rhs.data)[i];
			(*(data)).push_back(c);
			}
		}
	}

void Sequence::FromString(const string &Label, const string &Seq)
	{
	asserta(m_GSI == UINT_MAX);
	asserta(m_SMI == UINT_MAX);

	label = Label;
	data = new vector<char>;
	asserta(data);
	int L = (int) SIZE(Seq);
	(*(data)).push_back('@');
	for (int i = 0; i < L; ++i)
		{
		char c = Seq[i];
		(*(data)).push_back(c);
		}
	}

void Sequence::LogMe() const
	{
	Log("\n");
	uint L = GetLength();
	Log("Sequence(%p) length %u  >%s\n",
	  this, L, label.c_str());
	for (uint i = 1; i <= L; ++i)
		Log("%c", (*data)[i]);
	Log("\n");
	}
