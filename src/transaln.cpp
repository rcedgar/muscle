#include "muscle.h"
#include "transaln.h"

/***
PWPath:
	Pair-wise path Fresh + MSASeq {XYB}

MSAPath:
	Length = columns in MSA
	Path MSASeq to columns in MSA {MG}

TPath1:
	Fresh sequence F to MSA columns
	Length = columns in MSA plus nr. inserts in F {FGgI}
		F = fresh letter
		G = gap in F in pair-wise alignment
		I = insert in F in pair-wise alignment
		g = gap because MSA seq has gap in MSA

TPath2:
	Fresh sequence F to expanded MSA columns
	Length = columns in expanded MSA {FGgIi}
		F = fresh letter
		G = gap in F in pair-wise alignment
		g = gap because MSA seq has gap in MSA
		I = insert in F in pair-wise alignment
		i = padding insert because longer insert in other fresh sequence

MPath:
	MSA to expanded MSA
	Length = columns in expanded MSA {Mi}
		M = MSA column
		i = padding insert
***/

void TransAln::LogMe() const
	{
	Log("\n");
	Log("Pair-wise alignments:\n");
	const uint FreshCount = GetFreshCount();
	for (uint FreshIndex = 0; FreshIndex < FreshCount; ++FreshIndex)
		{
		const Sequence &FreshSeq = GetFreshSeq(FreshIndex);
		const uint MSAIndex = GetMSAIndex(FreshIndex);
		const Sequence &UngappedMSASeq = GetUngappedMSASeq(MSAIndex);
		const string &PWPath = GetPWPath(FreshIndex);
		LogAln(FreshSeq, UngappedMSASeq, PWPath);
		}

	Log("\n");
	Log("MSAPaths:\n");
	const uint MSACount = GetMSACount();
	for (uint MSAIndex = 0; MSAIndex < MSACount; ++MSAIndex)
		{
		const string &MSAPath = GetMSAPath(MSAIndex);
		const string &MSALabel = GetMSALabel(MSAIndex);
		Log("%s", MSAPath.c_str());
		Log("  >%s", MSALabel.c_str());
		Log("\n");
		}

	Log("\n");
	Log("MaxInserts:\n");
	uint M = 0;
	for (uint Col = 0; Col < SIZE(m_MSAColToMaxInserts); ++Col)
		{
		uint n = m_MSAColToMaxInserts[Col];
		if (n > 0)
			{
			Log(" [%u]=%u", Col, n);
			++M;
			}
		}
	Log(" (%u)\n", M);
	Log("ExtendedColCount = %u\n", m_ExtendedMSAColCount);
	Log("\n");

	Log("\n");
	Log("TPaths1:\n");
	for (uint i = 0; i < SIZE(m_TPaths1); ++i)
		{
		const string &Path1 = m_TPaths1[i];
		const string &FreshLabel = GetFreshLabel(i);
		Log("%s", Path1.c_str());
		Log("  >%s", FreshLabel.c_str());
		Log("\n");
		}

	for (uint i = 0; i < SIZE(m_TPaths1); ++i)
		LogTPath1Aln(i);

	Log("\n");
	Log("TPaths2:\n");
	for (uint i = 0; i < SIZE(m_TPaths2); ++i)
		{
		const string &Path2 = m_TPaths2[i];
		const string &FreshLabel = GetFreshLabel(i);
		Log("%s", Path2.c_str());
		Log("  >%s", FreshLabel.c_str());
		Log("\n");
		}
	Log("\n");
	Log("MPath\n");
	Log("%s\n", m_MPath.c_str());

	for (uint i = 0; i < SIZE(m_TPaths2); ++i)
		LogTPath2Aln(i, true);

	Log("\n");
	for (uint i = 0; i < MSACount; ++i)
		LogMPathAln(i, false);
	for (uint i = 0; i < SIZE(m_TPaths2); ++i)
		LogTPath2Aln(i, false);
	}

uint TransAln::GetFreshCount() const
	{
	asserta(m_FreshSeqs != 0);
	return m_FreshSeqs->GetSeqCount();
	}

uint TransAln::GetMSACount() const
	{
	asserta(m_MSA != 0);
	return m_MSA->GetSeqCount();
	}

uint TransAln::GetMSAIndex(uint FreshIndex) const
	{
	asserta(m_FreshIndexToMSAIndex != 0);
	asserta(FreshIndex < SIZE(*m_FreshIndexToMSAIndex));
	uint MSAIndex = (*m_FreshIndexToMSAIndex)[FreshIndex];
	return MSAIndex;
	}

const string &TransAln::GetPWPath(uint FreshIndex) const
	{
	asserta(m_PWPaths != 0);
	asserta(FreshIndex < SIZE(*m_PWPaths));
	const string &Path = (*m_PWPaths)[FreshIndex];
	return Path;
	}

const string &TransAln::GetMSAPath(uint FreshIndex) const
	{
	asserta(FreshIndex < SIZE(m_MSAPaths));
	return m_MSAPaths[FreshIndex];
	}

const string &TransAln::GetTPath1(uint FreshIndex) const
	{
	asserta(FreshIndex < SIZE(m_TPaths1));
	return m_TPaths1[FreshIndex];
	}

const string &TransAln::GetTPath2(uint FreshIndex) const
	{
	asserta(FreshIndex < SIZE(m_TPaths2));
	return m_TPaths2[FreshIndex];
	}

const string &TransAln::GetMSALabel(uint MSAIndex) const
	{
	const Sequence &Seq = GetMSASeq(MSAIndex);
	return Seq.GetLabel();
	}

const string &TransAln::GetFreshLabel(uint FreshIndex) const
	{
	const Sequence &Seq = GetFreshSeq(FreshIndex);
	return Seq.GetLabel();
	}

const Sequence &TransAln::GetMSASeq(uint MSAIndex) const
	{
	asserta(m_MSA != 0);
	const Sequence *Seq = m_MSA->GetSequence(MSAIndex);
	asserta(Seq != 0);
	return *Seq;
	}

const Sequence &TransAln::GetUngappedMSASeq(uint MSAIndex) const
	{
	asserta(MSAIndex < SIZE(m_UngappedMSASeqs));
	const Sequence *Seq = m_UngappedMSASeqs[MSAIndex];
	asserta(Seq != 0);
	return *Seq;
	}

const Sequence &TransAln::GetFreshSeq(uint FreshIndex) const
	{
	asserta(m_FreshSeqs != 0);
	const Sequence *Seq = m_FreshSeqs->GetSequence(FreshIndex);
	asserta(Seq != 0);
	return *Seq;
	}

uint TransAln::GetUngappedMSASeqLength(uint MSAIndex) const
	{
	const Sequence &Seq = GetUngappedMSASeq(MSAIndex);
	const uint L = Seq.GetLength();
	return L;
	}

uint TransAln::GetFreshSeqLength(uint FreshIndex) const
	{
	const Sequence &Seq = GetFreshSeq(FreshIndex);
	uint L = Seq.GetLength();
	return L;
	}

void TransAln::MakeTPath1(uint FreshIndex, string &Path1) const
	{
	Path1.clear();

	const uint MSAColCount = GetMSAColCount();
	const uint MSAIndex = GetMSAIndex(FreshIndex);
	const string &PWPath = GetPWPath(FreshIndex);
	const string &MSAPath = GetMSAPath(MSAIndex);
	const uint FreshL = GetFreshSeqLength(FreshIndex);
	const uint UL = GetUngappedMSASeqLength(MSAIndex);

	const uint PWColCount = SIZE(PWPath);
	uint NF = 0;
	uint NG = 0;
	uint NI = 0;
	uint Ng = 0;
	uint MSACol = 0;
	for (uint PWCol = 0; PWCol < PWColCount; ++PWCol)
		{
		char c = PWPath[PWCol];
		if (c == 'B' || c == 'Y')
			{
			while (MSAPath[MSACol] == 'G')
				{
				++MSACol;
				Path1 += 'g';
				++Ng;
				}
			}

		switch (c)
			{
		case 'B':
			Path1 += 'F';
			++NF;
			++MSACol;
			break;

		case 'X':
			Path1 += 'I';
			++NI;
			break;

		case 'Y':
			Path1 += 'G';
			++NG;
			++MSACol;
			break;

		default:
			asserta(false);
			}
		}

	while (MSACol < MSAColCount)
		{
		asserta(MSAPath[MSACol] == 'G');
		++MSACol;
		Path1 += 'g';
		++Ng;
		}

	asserta(NF + NG + Ng == MSAColCount);
	asserta(NF + NI == FreshL);
	asserta(NF + NG == UL);
	}

void TransAln::MakeMSAPath(uint MSAIndex, string &MSAPath) const
	{
	MSAPath.clear();
	const Sequence &MSASeq = GetMSASeq(MSAIndex);
	const byte *ByteSeq = MSASeq.GetBytePtr();
	asserta(MSASeq.GetLength() == m_MSAColCount);
	for (uint Col = 0; Col < m_MSAColCount; ++Col)
		{
		byte c = ByteSeq[Col];
		if (c == '-')
			MSAPath += 'G';
		else
			MSAPath += 'M';
		}
	}

void TransAln::MakeMPath(string &MPath) const
	{
	MPath.clear();
	asserta(SIZE(m_MSAColToMaxInserts) == m_MSAColCount + 1);
	for (uint Col = 0; Col <= m_MSAColCount; ++Col)
		{
		uint Ins = m_MSAColToMaxInserts[Col];
		for (uint i = 0; i < Ins; ++i)
			MPath += 'i';
		if (Col < m_MSAColCount)
			MPath += 'M';
		}
	}

void TransAln::Init(const MultiSequence &MSA,
  const MultiSequence &FreshSeqs,
  const vector<uint> &FreshIndexToMSAIndex,
  const vector<string> &PWPaths)
	{
	m_MSAPaths.clear();
	m_TPaths1.clear();
	m_TPaths2.clear();
	m_UngappedMSASeqs.clear();

	m_MSA = &MSA;
	m_FreshSeqs = &FreshSeqs;
	m_FreshIndexToMSAIndex = &FreshIndexToMSAIndex;
	m_PWPaths = &PWPaths;
	m_MSAColCount = MSA.GetColCount();

	const uint MSASeqCount = MSA.GetSeqCount();
	for (uint MSAIndex = 0; MSAIndex < MSASeqCount; ++MSAIndex)
		{
		string MSAPath;
		MakeMSAPath(MSAIndex, MSAPath);
		m_MSAPaths.push_back(MSAPath);

		const Sequence *Seq = MSA.GetSequence(MSAIndex);
		Sequence *UngappedSeq = Seq->CopyDeleteGaps();
		m_UngappedMSASeqs.push_back(UngappedSeq);
		}

	const uint FreshSeqCount = FreshSeqs.GetSeqCount();
	for (uint FreshIndex = 0; FreshIndex < FreshSeqCount; ++FreshIndex)
		{
		string Path1;
		MakeTPath1(FreshIndex, Path1);
		m_TPaths1.push_back(Path1);
		}

	SetMaxInserts();

	for (uint FreshIndex = 0; FreshIndex < FreshSeqCount; ++FreshIndex)
		{
		string Path2;
		MakeTPath2(FreshIndex, Path2);
		m_TPaths2.push_back(Path2);
		}
	MakeMPath(m_MPath);
	}

void TransAln::SetMaxInserts()
	{
	const uint FreshCount = GetFreshCount();
	asserta(SIZE(m_TPaths1) == FreshCount);
	m_MSAColToMaxInserts.clear();
	m_MSAColToMaxInserts.resize(m_MSAColCount+1, 0);
	for (uint FreshIndex = 0; FreshIndex < FreshCount; ++FreshIndex)
		{
		vector<uint> MSAColToInserts;
		MakeMSAColToInserts(FreshIndex, MSAColToInserts);
		asserta(SIZE(MSAColToInserts) == m_MSAColCount + 1);
		for (uint MSACol = 0; MSACol <= m_MSAColCount; ++MSACol)
			{
			uint Ins = MSAColToInserts[MSACol];
			m_MSAColToMaxInserts[MSACol] = max(Ins, m_MSAColToMaxInserts[MSACol]);
			}
		}

	m_ExtendedMSAColCount = 0;
	for (uint MSACol = 0; MSACol <= m_MSAColCount; ++MSACol)
		{
		uint Ins = m_MSAColToMaxInserts[MSACol];
		m_ExtendedMSAColCount += Ins;
		if (MSACol < m_MSAColCount)
			++m_ExtendedMSAColCount;
		}
	}

void TransAln::MakeMSAColToInserts(uint FreshIndex,
  vector<uint> &MSAColToInserts) const
	{
	MSAColToInserts.clear();
	const string &TPath1 = GetTPath1(FreshIndex);
	uint MSACol = 0;
	MSAColToInserts.resize(m_MSAColCount + 1, 0);
	const uint n = SIZE(TPath1);

	for (uint i = 0; i < n; ++i)
		{
		char c = TPath1[i];
		switch (c)
			{
		case 'F':
		case 'G':
		case 'g':
			++MSACol;
			break;

		case 'I':
			asserta(MSACol <= m_MSAColCount);
			MSAColToInserts[MSACol] += 1;
			break;

		default:
			asserta(false);
			}
		}

	asserta(MSACol == m_MSAColCount);
	}

void TransAln::MakeTPath2(uint FreshIndex, string &Path2) const
	{
	Path2.clear();
	asserta(SIZE(m_MSAColToMaxInserts) == m_MSAColCount + 1);
	const string &TPath1 = GetTPath1(FreshIndex);
	uint MSACol = 0;
	vector<uint> MSAColToInserts;
	MakeMSAColToInserts(FreshIndex, MSAColToInserts);
	const uint n1 = SIZE(TPath1);

	MSACol = 0;
	for (uint i = 0; i < n1; ++i)
		{
		char c = TPath1[i];
		Path2 += c;
		if (c != 'I')
			{
			asserta(MSACol < m_MSAColCount);
			uint InsertCount = MSAColToInserts[MSACol];
			uint MaxInsertCount = m_MSAColToMaxInserts[MSACol];
			asserta(InsertCount <= MaxInsertCount);
			for (uint j = InsertCount; j < MaxInsertCount; ++j)
				Path2 += 'i';
			}

		switch (c)
			{
		case 'F':
		case 'G':
		case 'g':
			++MSACol;
			break;

		case 'I':
			break;

		default:
			asserta(false);
			}
		}
	asserta(MSACol == m_MSAColCount);
	uint InsertCount = MSAColToInserts[m_MSAColCount];
	uint MaxInsertCount = m_MSAColToMaxInserts[m_MSAColCount];
	asserta(InsertCount <= MaxInsertCount);
	for (uint j = InsertCount; j < MaxInsertCount; ++j)
		Path2 += 'i';

	if (SIZE(Path2) != m_ExtendedMSAColCount)
		{
		LogMe();
		Log("FreshIndex %u, MSAIndex %u, Path2=%s\n",
		  FreshIndex, GetMSAIndex(FreshIndex), Path2.c_str());
		Die("|Path2|=%u, m_ExtendedMSAColCount=%u",
		  SIZE(Path2), m_ExtendedMSAColCount);
		}
	}

void TransAln::LogTPath1Aln(uint FreshIndex) const
	{
	const string &TPath1 = GetTPath1(FreshIndex);
	const uint MSAIndex = GetMSAIndex(FreshIndex);
	const Sequence &F = GetFreshSeq(FreshIndex);
	const Sequence &U = GetUngappedMSASeq(MSAIndex);
	const uint FL = F.GetLength();
	const uint UL = U.GetLength();
	const string &FLabel = GetFreshLabel(FreshIndex);
	const string &ULabel = GetMSALabel(MSAIndex);
	const byte *Fb = F.GetBytePtr();
	const byte *Ub = U.GetBytePtr();
	const uint ColCount = SIZE(TPath1);
	uint FPos = 0;
	uint UPos = 0;
	string FRow;
	string URow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = TPath1[Col];
		switch (c)
			{
		case 'F':
			{
			char f = Fb[FPos];
			char u = Ub[UPos];
			FRow += f;
			URow += u;
			++FPos;
			++UPos;
			break;
			}

		case 'G':
			{
			char u = Ub[UPos];
			FRow += '-';
			URow += u;
			++UPos;
			break;
			}

		case 'I':
			{
			char f = Fb[FPos];
			FRow += f;
			URow += '.';
			++FPos;
			break;
			}

		case 'g':
			{
			FRow += '.';
			URow += '.';
			break;
			}
			}
		}
	Log("\n");
	Log("%s\n", TPath1.c_str());
	Log("%s  >%s\n", FRow.c_str(), FLabel.c_str());
	Log("%s  >%s\n", URow.c_str(), ULabel.c_str());
	asserta(FPos == FL);
	asserta(UPos == UL);
	}

void TransAln::LogTPath2Aln(uint FreshIndex, bool WithPath) const
	{
	const string &TPath2 = GetTPath2(FreshIndex);
	const uint MSAIndex = GetMSAIndex(FreshIndex);
	const Sequence &F = GetFreshSeq(FreshIndex);
	const Sequence &U = GetUngappedMSASeq(MSAIndex);
	const uint FL = F.GetLength();
	const uint UL = U.GetLength();
	const string &FLabel = GetFreshLabel(FreshIndex);
	const string &ULabel = GetMSALabel(MSAIndex);
	const byte *Fb = F.GetBytePtr();
	const byte *Ub = U.GetBytePtr();
	const uint ColCount = SIZE(TPath2);
	uint FPos = 0;
	uint UPos = 0;
	string FRow;
	string URow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = TPath2[Col];
		switch (c)
			{
		case 'F':
			{
			char f = Fb[FPos];
			char u = Ub[UPos];
			FRow += f;
			URow += u;
			++FPos;
			++UPos;
			break;
			}

		case 'G':
			{
			char u = Ub[UPos];
			FRow += '-';
			URow += u;
			++UPos;
			break;
			}

		case 'I':
			{
			char f = Fb[FPos];
			FRow += f;
			URow += '.';
			++FPos;
			break;
			}

		case 'g':
			{
			FRow += '.';
			URow += '.';
			break;
			}

		case 'i':
			{
			FRow += '.';
			URow += '.';
			break;
			}
			}
		}
	if (WithPath)
		{
		Log("\n");
		Log("%s\n", TPath2.c_str());
		}
	Log("%s  [F] >%s\n", FRow.c_str(), FLabel.c_str());
	Log("%s  [U] >%s\n", URow.c_str(), ULabel.c_str());
	asserta(FPos == FL);
	asserta(UPos == UL);
	}

void TransAln::LogMPathAln(uint MSAIndex, bool WithPath) const
	{
	const Sequence &M = GetMSASeq(MSAIndex);
	const byte *Mb = M.GetBytePtr();
	const uint ColCount = SIZE(m_MPath);
	const string &MLabel = GetMSALabel(MSAIndex);
	uint MSACol = 0;
	string MRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_MPath[Col];
		switch (c)
			{
		case 'M':
			{
			MRow += Mb[MSACol];
			++MSACol;
			break;
			}

		case 'i':
			{
			MRow += '.';
			break;
			}
			}
		}
	if (WithPath)
		{
		Log("\n");
		Log("%s\n", m_MPath.c_str());
		}
	Log("%s  [M] >%s\n", MRow.c_str(), MLabel.c_str());
	asserta(MSACol == m_MSAColCount);
	}

Sequence *TransAln::ExtendFreshSeq(uint FreshIndex) const
	{
	const Sequence &F = GetFreshSeq(FreshIndex);
	const byte *Fb = F.GetBytePtr();
	const string &TPath2 = GetTPath2(FreshIndex);
	const uint ColCount = SIZE(m_MPath);
	const string &FLabel = GetFreshLabel(FreshIndex);
	uint MSACol = 0;
	Sequence *FX = Sequence::NewSequence();
	FX->InitData();
	FX->m_Label = F.m_Label;
	uint FPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = TPath2[Col];
		switch (c)
			{
		case 'F':
		case 'I':
			{
			char f = Fb[FPos];
			FX->AppendChar(f);
			++FPos;
			break;
			}

		case 'G':
		case 'g':
		case 'i':
			{
			FX->AppendChar('-');
			break;
			}
		
		default:
			Die("Invalid char '%c' in TPath2", c);
			}
		}
	asserta(FX->GetLength() == m_ExtendedMSAColCount);
	return FX;
	}

Sequence *TransAln::ExtendMSASeq(uint MSAIndex) const
	{
	const Sequence &M = GetMSASeq(MSAIndex);
	const byte *Mb = M.GetBytePtr();
	const uint ColCount = SIZE(m_MPath);
	const string &MLabel = GetMSALabel(MSAIndex);
	uint MSACol = 0;
	Sequence *MX = Sequence::NewSequence();
	MX->InitData();
	MX->m_Label = M.m_Label;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = m_MPath[Col];
		switch (c)
			{
		case 'M':
			{
			char c = Mb[MSACol];
			MX->AppendChar(c);
			++MSACol;
			break;
			}

		case 'i':
			{
			MX->AppendChar('-');
			break;
			}

		default:
			Die("Invalid char '%c' in MPath", c);
			}
		}
	asserta(MX->GetLength() == m_ExtendedMSAColCount);
	return MX;
	}

void TransAln::MakeExtendedMSA()
	{
	const uint MSACount = GetMSACount();
	const uint FreshCount = GetFreshCount();
	m_ExtendedMSA = new MultiSequence;

	for (uint i = 0; i < MSACount; ++i)
		{
		Sequence *S = ExtendMSASeq(i);
		m_ExtendedMSA->AddSequence(S, true);
		}

	for (uint i = 0; i < FreshCount; ++i)
		{
		Sequence *S = ExtendFreshSeq(i);
		m_ExtendedMSA->AddSequence(S, true);
		}
	}

void cmd_transaln()
	{
	const string &InputFileName = opt(transaln);
	const string &RefFileName = opt(ref);
	const string &OutputFileName = opt(output);

	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName);
	const uint InputSeqCount = InputSeqs.GetSeqCount();

	MultiSequence RefMSA;
	RefMSA.FromFASTA(RefFileName);

	const uint RefSeqCount = RefMSA.GetSeqCount();
	MultiSequence UngappedRefSeqs;
	for (uint i = 0; i < RefSeqCount; ++i)
		{
		const Sequence *AlignedRefSeq = RefMSA.GetSequence(i);
		Sequence *UngappedRefSeq = AlignedRefSeq->CopyDeleteGaps();
		UngappedRefSeqs.AddSequence(UngappedRefSeq, true);
		}

	vector<string> PWPaths;
	vector<uint> FreshIndexToMSAIndex;
	for (uint InputSeqIndex = 0; InputSeqIndex < InputSeqCount; ++InputSeqIndex)
		{
		const uint RefSeqIndex = InputSeqIndex%RefSeqCount;
		FreshIndexToMSAIndex.push_back(RefSeqIndex);

		const Sequence *InputSeq = InputSeqs.GetSequence(InputSeqIndex);
		const Sequence *RefSeq = UngappedRefSeqs.GetSequence(RefSeqIndex);
		const string &InputLabel = InputSeq->m_Label;
		const string &RefLabel = RefSeq->m_Label;

		string PWPath;
		AlignPairFlat(InputLabel, RefLabel, PWPath);

		PWPaths.push_back(PWPath);
		}
	TransAln TA;
	TA.Init(RefMSA, InputSeqs, FreshIndexToMSAIndex, PWPaths);
	TA.MakeExtendedMSA();
	asserta(TA.m_ExtendedMSA != 0);
	TA.m_ExtendedMSA->WriteMFA(OutputFileName);
	}
