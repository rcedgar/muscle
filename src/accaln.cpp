#include "muscle.h"
#include "accaln.h"

/***

0       1j46_A  1aab_   0.498   15.6%
0       MQ---DRVKRPMNAFIVWSRDQRRKMALENPR--MRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK
0       GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTY-------IPPKGE
0       21   112223333333333333333333332  3445678889999999999999998776655443322222111       112236
0       27   472680898876666667542111007  6065750460368999988776317517371737073210875       360882

1       1j46_A  1k99_A  0.566   19.8%
1       MQ------DRVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK
1       MKKLKKHPDFPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK
1       31      22455566655555555555555555555667888999999999999999987766655443322222211111111222346
1       19      38116910099888889987654443567066269035799999887774298496272713075432094322358269559

***/

void AssertSameSeq(const string &Row1, const string &Row2)
	{
	string Seq1;
	string Seq2;

	for (uint i = 0; i < SIZE(Row1); ++i)
		{
		char c = Row1[i];
		if (!isgap(c))
			Seq1 += toupper(c);
		}

	for (uint i = 0; i < SIZE(Row2); ++i)
		{
		char c = Row2[i];
		if (!isgap(c))
			Seq2 += toupper(c);
		}
	asserta(Seq1 == Seq2);
	}

void AccAln::FromFile(const string &FileName)
	{
	Clear();
	m_FileName = FileName;
	m_f = OpenStdioFile(m_FileName);
	m_PairIndex = 0;
	for (;;)
		{
		bool Ok = NextPair();
		if (!Ok)
			break;
		++m_PairIndex;
		}

	CloseStdioFile(m_f);
	m_f = 0;
	}

bool AccAln::GetBlankLine()
	{
	++m_LineNr;
	string Line;
	bool Ok = ReadLineStdioFile(m_f, Line);
	if (!Ok)
		return false;
	StripWhiteSpace(Line);
	if (!Line.empty())
		Die("%s:%u expected blank line, got '%s'",
		  m_FileName.c_str(), m_LineNr, Line.c_str());
	return true;
	}

void AccAln::GetAlnRow(string &Row)
	{
	Row.clear();

	++m_LineNr;
	string Line;
	bool Ok = ReadLineStdioFile(m_f, Line);
	if (!Ok)
		Die("%s:%u unexpected end-of-file in pair %u",
		  m_FileName.c_str(), m_LineNr, m_PairIndex);

	const char *ptrLine = Line.c_str();
	const char *ptrTab = strchr(ptrLine, '\t');
	if (ptrTab == 0)
		Die("%s:%u no tab in '%s'",
		  m_FileName.c_str(), m_LineNr, Line.c_str());

	string StrInt;
	for (const char *p = ptrLine; p != ptrTab; ++p)
		StrInt += *p;
	uint PairIndex = StrToUint(StrInt);
	if (PairIndex != m_PairIndex)
		Die("%s:%u expected pair %u got %u in '%s'",
		  m_FileName.c_str(), m_LineNr, m_PairIndex, PairIndex, Line.c_str());

	for (const char *p = ptrTab + 1; *p; ++p)
		Row += *p;
	}

bool AccAln::NextPair()
	{
	string BlankLine;
	bool Ok = GetBlankLine();
	if (!Ok)
		return false;

	string Row1;
	string Row2;
	string Row3;
	string Row4;
	string Row5;
	GetAlnRow(Row1);
	GetAlnRow(Row2);
	GetAlnRow(Row3);
	GetAlnRow(Row4);
	GetAlnRow(Row5);

	const uint ColCount = SIZE(Row2);
	if (SIZE(Row3) != ColCount ||
		SIZE(Row4) != ColCount ||
		SIZE(Row5) != ColCount)
			Die("Rows different lengths pair %u", m_PairIndex);

	vector<string> Fields;
	Split(Row1, Fields, '\t');
	if (SIZE(Fields) != 4)
		Die("Expected 4 fields in '%s'", Row1.c_str());
	
	const string &Label1 = Fields[0];
	const string &Label2 = Fields[1];
	double EA = StrToFloat(Fields[2]);

	vector<float> EAs;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char d1 = Row4[Col];
		char d2 = Row5[Col];
		float EA = -1;
		if (d1 == ' ' && d2 == ' ')
			{
			char c1 = Row2[Col];
			char c2 = Row3[Col];
			if (!(isgap(c1) || isgap(c2)))
				Die("%s:%u col %u expected gap",
				  m_FileName.c_str(), m_LineNr, Col);
			EA = -1;
			}
		else if (d1 == '^' && d2 == '^')
			EA = 1;
		else if (isdigit(d1) && isdigit(d2))
			{
			char Str[3];
			Str[0] = d1;
			Str[1] = d2;
			Str[2] = 0;
			int i = atoi(Str);
			asserta(i >= 0 && i <= 99);
			EA = float(i)/100;
			}
		else
			Die("%s:%u col %u invalid EA %c%c", 
			  m_FileName.c_str(), m_LineNr, Col, d1, d2);
		EAs.push_back(EA);
		}

	m_Labels1.push_back(Label1);
	m_Labels2.push_back(Label2);
	m_Rows1.push_back(Row2);
	m_Rows2.push_back(Row3);
	m_EADigits1.push_back(Row4);
	m_EADigits2.push_back(Row5);
	m_EAs.push_back(EAs);
	return true;
	}

uint AccAln::GetPairCount() const
	{
	uint PairCount = SIZE(m_Labels1);
	asserta(SIZE(m_Labels2) == PairCount);
	asserta(SIZE(m_Rows1) == PairCount);
	asserta(SIZE(m_Rows2) == PairCount);
	asserta(SIZE(m_EAs) == PairCount);
	return PairCount;
	}

void AccAln::ReadRef(const string &FileName)
	{
	m_RefFileName = FileName;
	m_RefMSA.FromFASTAFile_PreserveCase(FileName);
	const uint RefSeqCount = m_RefMSA.GetSeqCount();

	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		string Label = m_RefMSA.GetSeqName(RefSeqIndex);
		if (m_RefLabelToIndex.find(Label) != m_RefLabelToIndex.end())
			Die("Dupe ref label >%s", Label.c_str());
		m_RefLabelToIndex[Label] = RefSeqIndex;
		}
	}

void AccAln::AlnRefPair(uint PairIndex)
	{
	const string &Label1 = m_Labels1[PairIndex];
	const string &Label2 = m_Labels2[PairIndex];
	if (m_RefLabelToIndex.find(Label1) == m_RefLabelToIndex.end())
		Die("Test label not found >%s", Label1.c_str());
	if (m_RefLabelToIndex.find(Label2) == m_RefLabelToIndex.end())
		Die("Test label not found >%s", Label2.c_str());
	const uint RefSeqIndex1 = m_RefLabelToIndex[Label1];
	const uint RefSeqIndex2 = m_RefLabelToIndex[Label2];

	string RefRow1;
	string RefRow2;
	m_RefMSA.GetRowStr(RefSeqIndex1, RefRow1);
	m_RefMSA.GetRowStr(RefSeqIndex2, RefRow2);

	const string &TestRow1 = m_Rows1[PairIndex];
	const string &TestRow2 = m_Rows2[PairIndex];
	const vector<float> &EAs = m_EAs[PairIndex];

	AssertSameSeq(TestRow1, RefRow1);
	AssertSameSeq(TestRow2, RefRow2);

	m_RefRows1.push_back(RefRow1);
	m_RefRows2.push_back(RefRow2);
	}

void AccAln::GetCorrect(uint PairIndex, vector<bool> &ColToCorrect,
  vector<bool> &ColToAligned)
	{
	ColToCorrect.clear();
	ColToAligned.clear();

	const string &RefRow1 = m_RefRows1[PairIndex];
	const string &RefRow2 = m_RefRows2[PairIndex];
	const uint RefColCount = SIZE(RefRow1);
	asserta(SIZE(RefRow2) == RefColCount);

	vector<uint> PosToRefCol1;
	vector<uint> PosToRefCol2;
	vector<bool> PosToAligned1;
	vector<bool> PosToAligned2;

	for (uint Col = 0; Col < RefColCount; ++Col)
		{
		char c1 = RefRow1[Col];
		char c2 = RefRow2[Col];
		bool Gap1 = isgap(c1);
		bool Gap2 = isgap(c2);
		if (Gap1 && Gap2)
			continue;

		bool Aligned = (isupper(c1) && isupper(c2));
		if (!Gap1)
			{
			PosToRefCol1.push_back(Col);
			PosToAligned1.push_back(Aligned);
			}
		if (!Gap2)
			{
			PosToRefCol2.push_back(Col);
			PosToAligned2.push_back(Aligned);
			}
		}

	const uint L1 = SIZE(PosToRefCol1);
	const uint L2 = SIZE(PosToRefCol2);

	const string &TestRow1 = m_Rows1[PairIndex];
	const string &TestRow2 = m_Rows2[PairIndex];
	const uint TestColCount = SIZE(TestRow1);
	asserta(SIZE(TestRow2) == TestColCount);
	uint Pos1 = 0;
	uint Pos2 = 0;
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		char c1 = TestRow1[Col];
		char c2 = TestRow2[Col];
		bool Gap1 = isgap(c1);
		bool Gap2 = isgap(c2);
		
		if (Gap1 && Gap2)
			continue;

		if (Gap1 || Gap2)
			{
			ColToCorrect.push_back(false);
			ColToAligned.push_back(false);
			}
		else
			{
			asserta(isalpha(c1) && isalpha(c2));
			asserta(Pos1 < L1);
			asserta(Pos2 < L2);
			uint RefCol1 = PosToRefCol1[Pos1];
			uint RefCol2 = PosToRefCol2[Pos2];
			bool Correct = (RefCol1 == RefCol2);
			asserta(Pos1 < SIZE(PosToAligned1));
			asserta(Pos2 < SIZE(PosToAligned2));
			bool Aligned = (PosToAligned1[Pos1] || PosToAligned2[Pos2]);
			ColToAligned.push_back(Aligned);
			ColToCorrect.push_back(Correct);
			}

		if (!Gap1)
			++Pos1;
		if (!Gap2)
			++Pos2;
		}
	}

void AccAln::AlnRef()
	{
	const uint PairCount = GetPairCount();
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		AlnRefPair(PairIndex);

	m_TestColToCorrectVec.clear();
	m_TestColToAlignedVec.clear();
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		vector<bool> ColToCorrect;
		vector<bool> ColToAligned;
		GetCorrect(PairIndex, ColToCorrect, ColToAligned);
		m_TestColToCorrectVec.push_back(ColToCorrect);
		m_TestColToAlignedVec.push_back(ColToAligned);
		}
	}

float AccAln::GetAlignedEA(uint PairIndex) const
	{
	const vector<float> &EAs = m_EAs[PairIndex];
	const vector<bool> &TestColToAligned = m_TestColToAlignedVec[PairIndex];
	const uint TestColCount = SIZE(EAs);
	asserta(SIZE(TestColToAligned) == TestColCount);

	uint N = 0;
	float SumEA = 0;
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		if (TestColToAligned[Col])
			{
			++N;
			SumEA += EAs[Col];
			}
		}
	if (N == 0)
		return 0;
	float AvgEA = SumEA/N;
	return AvgEA;
	}

float AccAln::GetQ(uint PairIndex) const
	{
	const vector<float> &EAs = m_EAs[PairIndex];
	const vector<bool> &TestColToAligned = m_TestColToAlignedVec[PairIndex];
	const vector<bool> &TestColToCorrect = m_TestColToCorrectVec[PairIndex];
	const uint TestColCount = SIZE(EAs);
	asserta(SIZE(TestColToAligned) == TestColCount);
	asserta(SIZE(TestColToCorrect) == TestColCount);

	uint N = 0;
	uint n = 0;
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		if (TestColToAligned[Col])
			{
			++N;
			if (TestColToCorrect[Col])
				++n;
			}
		}
	if (N == 0)
		return 0;
	float Q = float(n)/N;
	return Q;
	}

void AccAln::WriteRefAln(FILE *f, uint PairIndex) const
	{
	if (f == 0)
		return;

	const string &Label1 = m_Labels1[PairIndex];
	const string &Label2 = m_Labels2[PairIndex];
	const string &TestRow1 = m_Rows1[PairIndex];
	const string &TestRow2 = m_Rows2[PairIndex];
	const string &EADigits1 = m_EADigits1[PairIndex];
	const string &EADigits2 = m_EADigits2[PairIndex];
	const vector<float> &EAs = m_EAs[PairIndex];
	const vector<bool> &TestColToAligned = m_TestColToAlignedVec[PairIndex];
	const vector<bool> &TestColToCorrect = m_TestColToCorrectVec[PairIndex];
	const uint TestColCount = SIZE(TestRow1);
	asserta(SIZE(TestRow2) == TestColCount);
	asserta(SIZE(EAs) == TestColCount);
	asserta(SIZE(TestColToAligned) == TestColCount);
	asserta(SIZE(TestColToCorrect) == TestColCount);
	asserta(SIZE(EADigits1) == TestColCount);
	asserta(SIZE(EADigits2) == TestColCount);

	const string &RefRow1 = m_RefRows1[PairIndex];
	const string &RefRow2 = m_RefRows2[PairIndex];
	const uint RefColCount = SIZE(RefRow1);
	asserta(SIZE(RefRow2) == RefColCount);

	float EA = GetAlignedEA(PairIndex);
	float Q = GetQ(PairIndex);

	fprintf(f, "\n");
	fprintf(f, "%u.info	%s	%s	EA=%.3g	Q=%.3g\n",
	  PairIndex, Label1.c_str(), Label2.c_str(), EA, Q);

	fprintf(f, "%u.ref1	", PairIndex);
	for (uint Col = 0; Col < RefColCount; ++Col)
		{
		char c1 = RefRow1[Col];
		char c2 = RefRow2[Col];
		if (isgap(c1) && isgap(c2))
			continue;
		fputc(c1, f);
		}
	fprintf(f, "\n");

	fprintf(f, "%u.ref2	", PairIndex);
	for (uint Col = 0; Col < RefColCount; ++Col)
		{
		char c1 = RefRow1[Col];
		char c2 = RefRow2[Col];
		if (isgap(c1) && isgap(c2))
			continue;
		fputc(c2, f);
		}
	fprintf(f, "\n");
	fprintf(f, "\n");

	fprintf(f, "%u.tst1	", PairIndex);
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		char c1 = TestRow1[Col];
		char c2 = TestRow2[Col];
		if (isgap(c1) && isgap(c2))
			continue;
		bool Aligned = TestColToAligned[Col];
		if (Aligned)
			c1 = toupper(c1);
		else
			c1 = tolower(c1);
		fputc(c1, f);
		}
	fprintf(f, "\n");

	fprintf(f, "%u.tst2	", PairIndex);
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		char c1 = TestRow1[Col];
		char c2 = TestRow2[Col];
		if (isgap(c1) && isgap(c2))
			continue;
		bool Aligned = TestColToAligned[Col];
		if (Aligned)
			c2 = toupper(c2);
		else
			c2 = tolower(c2);
		fputc(c2, f);
		}
	fprintf(f, "\n");

	fprintf(f, "%u.ead1	%s\n", PairIndex, EADigits1.c_str());
	fprintf(f, "%u.ead2	%s\n", PairIndex, EADigits2.c_str());

	fprintf(f, "%u.corr	", PairIndex);
	for (uint Col = 0; Col < TestColCount; ++Col)
		{
		char c1 = TestRow1[Col];
		char c2 = TestRow2[Col];
		if (isgap(c1) && isgap(c2))
			continue;
		bool Aligned = TestColToAligned[Col];
		bool Correct = TestColToCorrect[Col];
		if (!Aligned)
			fputc('.', f);
		else
			fputc(Correct ? '|' : 'x', f);
		}
	fprintf(f, "\n");
	}

void AccAln::WriteRefAln(const string &FileName) const
	{
	if (FileName.empty())
		return;

	FILE *f = CreateStdioFile(FileName);
	const uint PairCount = GetPairCount();
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		WriteRefAln(f, PairIndex);
	CloseStdioFile(f);
	}

void cmd_accalnref()
	{
	const string &AAFileName = opt(accalnref);
	const string &RefFileName = opt(ref);
	const string &OutputFileName = opt(output);

	AccAln AA;
	AA.FromFile(AAFileName);
	AA.ReadRef(RefFileName);
	AA.AlnRef();
	AA.WriteRefAln(OutputFileName);
	}
