#include "muscle.h"
#include "pprog.h"

void ReadStringsFromFile(const string &FileName,
  vector<string> &Strings)
	{
	Strings.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		Strings.push_back(Line);
	}

void InvertPath(const string &Path, string &InvertedPath)
	{
	InvertedPath.clear();
	const uint n = SIZE(Path);
	for (uint i = 0; i < n; ++i)
		{
		char c = Path[i];
		if (c == 'B')
			InvertedPath.push_back('B');
		else if (c == 'X')
			InvertedPath.push_back('Y');
		else if (c == 'Y')
			InvertedPath.push_back('X');
		else
			Die("Invalid path char '%c'", c);
		}
	}

void ValidatePath(const string &Path, uint LX, uint LY)
	{
	uint nX = 0;
	uint nY = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = toupper(Path[i]);
		switch (c)
			{
		case 'X':	++nX; break;
		case 'Y':	++nY; break;
		case 'B':	++nX, ++nY; break;
		default:	asserta(false);
			}
		}
	asserta(nX == LX);
	asserta(nY == LY);
	}

void AlignMSAsByPath(const MultiSequence &MSA1, const MultiSequence &MSA2,
  const string &Path, MultiSequence &MSA12)
	{
	uint LX = MSA1.GetColCount();
	uint LY = MSA2.GetColCount();
	ValidatePath(Path, LX, LY);
	const uint SeqCount1 = MSA1.GetSeqCount();
	const uint SeqCount2 = MSA2.GetSeqCount();

	for (int SeqIndex = 0; SeqIndex < (int) SeqCount1; ++SeqIndex)
		{
		const Sequence *Seq1 = MSA1.GetSequence(SeqIndex);
		Sequence *AlignedSeq1 = Seq1->AddGapsPath(Path, 'X');
		MSA12.AddSequence(AlignedSeq1, true);
		}

	for (int SeqIndex = 0; SeqIndex < MSA2.GetNumSequences(); ++SeqIndex)
		{
		const Sequence *Seq2 = MSA2.GetSequence(SeqIndex);
		Sequence *AlignedSeq2 = Seq2->AddGapsPath(Path, 'Y');
		MSA12.AddSequence(AlignedSeq2, true);
		}
	}

void PProg::DeleteIndexesFromPending(uint Index1, uint Index2)
	{
	bool Found1 = false;
	bool Found2 = false;
	vector<uint> NewPending;
	for (uint i = 0; i < SIZE(m_Pending); ++i)
		{
		uint Index = m_Pending[i];
		if (Index == Index1)
			{
			asserta(!Found1);
			Found1 = true;
			continue;
			}
		if (Index == Index2)
			{
			asserta(!Found2);
			Found2 = true;
			continue;
			}
		NewPending.push_back(Index);
		}
	asserta(Found1);
	asserta(Found2);
	asserta(SIZE(NewPending) + 2 == SIZE(m_Pending));
	m_Pending = NewPending;
	}

void PProg::FindBestPair(uint &BestIndex1, uint &BestIndex2) const
	{
	const uint N = SIZE(m_Pending);
	asserta(N >= 2);
	BestIndex1 = m_Pending[0];
	BestIndex2 = m_Pending[1];
	asserta(BestIndex1 < m_NodeCount);
	asserta(BestIndex2 < m_NodeCount);
	float BestScore = m_ScoreMx[BestIndex1][BestIndex2];
	for (uint i = 0; i < N; ++i)
		{
		uint Indexi = m_Pending[i];
		for (uint j = i+1; j < N; ++j)
			{
			uint Indexj = m_Pending[j];
			asserta(Indexi < SIZE(m_ScoreMx));
			asserta(Indexj < SIZE(m_ScoreMx[Indexi]));
			float Score = m_ScoreMx[Indexi][Indexj];
			if (Score > BestScore)
				{
				BestScore = Score;
				BestIndex1 = Indexi;
				BestIndex2 = Indexj;
				}
			}
		}
	}

const string &PProg::GetMSALabel(uint Index) const
	{
	asserta(Index < SIZE(m_MSALabels));
	const string &Label = m_MSALabels[Index];
	asserta(!Label.empty());
	return Label;
	}

void PProg::SetMSA(uint Index, const MultiSequence &MSA)
	{
	asserta(Index < SIZE(m_MSAs));
//	asserta(m_MSAs[Index] == 0); // TODO: memory leak
	m_MSAs[Index] = &MSA;
	}

void PProg::SetMSALabel(uint Index, const string &Label)
	{
	asserta(Index < SIZE(m_MSALabels));
	asserta(m_MSALabels[Index].empty());
	m_MSALabels[Index] = Label;
	asserta(m_MSALabelToIndex.find(Label) == m_MSALabelToIndex.end());
	m_MSALabelToIndex[Label] = Index;
	}

const MultiSequence &PProg::GetMSA(uint Index) const
	{
	asserta(Index < SIZE(m_MSAs));
	const MultiSequence *MSA = m_MSAs[Index];
	asserta(MSA != 0);
	return *MSA;
	}

const MultiSequence &PProg::GetFinalMSA() const
	{
	asserta(m_InputMSACount > 0);
	uint FinalIndex = 2*(m_InputMSACount - 1);
	asserta(FinalIndex + 1 == SIZE(m_MSAs));
	const MultiSequence *MSA = m_MSAs[FinalIndex];
	asserta(MSA != 0);
	return *MSA;
	}

void PProg::SetMSAs(const vector<const MultiSequence *> &MSAs,
  const vector<string> &MSALabels)
	{
	m_MSALabelToIndex.clear();
	m_InputMSACount = SIZE(MSAs);
	asserta(SIZE(MSALabels) == m_InputMSACount);
	m_MSAs = MSAs;
	m_MSALabels = MSALabels;

	uint TotalMSACount = 2*m_InputMSACount - 1;
	m_MSAs.resize(TotalMSACount, 0);
	m_MSALabels.resize(TotalMSACount, "");

	for (uint MSAIndex = 0; MSAIndex < m_InputMSACount; ++MSAIndex)
		{
		const string &MSALabel = MSALabels[MSAIndex];
		asserta(m_MSALabelToIndex.find(MSALabel) == m_MSALabelToIndex.end());
		m_MSALabelToIndex[MSALabel] = MSAIndex;
		}
	}

void PProg::LoadMSAs(const vector<string> &FileNames, bool &IsNucleo)
	{
	m_MSALabelToIndex.clear();
	m_MSAs.clear();
	m_MSALabels.clear();

	m_InputMSACount = SIZE(FileNames);
	asserta(m_InputMSACount > 1);

	uint TotalMSACount = 2*m_InputMSACount - 1;
	m_MSAs.resize(TotalMSACount, 0);
	m_MSALabels.resize(TotalMSACount, "");

	m_JoinCount = m_InputMSACount - 1;
	m_NodeCount = m_InputMSACount + m_JoinCount;

	for (uint MSAIndex = 0; MSAIndex < m_InputMSACount; ++MSAIndex)
		{
		const string &FileName = FileNames[MSAIndex];
		ProgressStep(MSAIndex, m_InputMSACount, "Reading %s", FileName.c_str());
		MultiSequence &MSA = *new MultiSequence;
		MSA.LoadMFA(FileName);
		bool IsNuc = MSA.GuessIsNucleo();
		if (MSAIndex == 0)
			IsNucleo = IsNuc;
		else
			asserta(IsNucleo == IsNuc);

		string MSALabel;
		GetBaseName(FileName, MSALabel);

		SetMSALabel(MSAIndex, MSALabel);
		SetMSA(MSAIndex, MSA);
		}
	}

void PProg::Run()
	{
	m_JoinMSAIndexes1.clear();
	m_JoinMSAIndexes2.clear();
	m_ScoreMx.clear();
	m_PathMx.clear();

	for (uint i = 0; i < m_InputMSACount; ++i)
		m_Pending.push_back(i);

	m_JoinCount = m_InputMSACount - 1;
	m_NodeCount = m_InputMSACount + m_JoinCount;
	AlignAllInputPairs();

	for (m_JoinIndex = 0; m_JoinIndex < m_JoinCount; ++m_JoinIndex)
		{
		ProgressLog("____________________________________________\n");
		ProgressLog("Join %u/%u, pending %u\n", 
		  m_JoinIndex+1, m_JoinCount, SIZE(m_Pending));
		uint Index1;
		uint Index2;
		FindBestPair(Index1, Index2);
		asserta(Index1 != Index2);
		Join_ByPrecomputedPath(Index1, Index2);
		AlignNewToPending();
		}
	}

void PProg::AlignAllInputPairs()
	{
	const uint PairCount = (m_InputMSACount*(m_InputMSACount - 1))/2;
	m_ScoreMx.resize(m_NodeCount);
	m_PathMx.resize(m_NodeCount);
	for (uint i = 0; i < m_NodeCount; ++i)
		{
		m_ScoreMx[i].resize(m_NodeCount);
		m_PathMx[i].resize(m_NodeCount);
		}

	uint PairIndex = 0;
	for (uint MSAIndex1 = 0; MSAIndex1 < m_InputMSACount; ++MSAIndex1)
		{
		const string &MSALabel1 = GetMSALabel(MSAIndex1);
		const MultiSequence &MSA1 = GetMSA(MSAIndex1);
		const uint SeqCount1 = MSA1.GetSeqCount();
		for (uint MSAIndex2 = MSAIndex1+1; MSAIndex2 < m_InputMSACount; ++MSAIndex2)
			{
			++PairIndex;
			Progress("Input pair %u / %u (%.1f%%)\n",
			  PairIndex, PairCount, GetPct(PairIndex, PairCount));

			const string &MSALabel2 = GetMSALabel(MSAIndex2);
			const MultiSequence &MSA2 = GetMSA(MSAIndex2);
			const uint SeqCount2 = MSA2.GetSeqCount();

			string Path;
			float Score = AlignMSAsFlat(MSALabel1 + "+" + MSALabel2,
			  MSA1, MSA2, m_TargetPairCount, Path);

			const uint ColCount1 = MSA1.GetColCount();
			const uint ColCount2 = MSA2.GetColCount();
			ValidatePath(Path, ColCount1, ColCount2);

			string InvertedPath;
			InvertPath(Path, InvertedPath);
			ValidatePath(InvertedPath, ColCount2, ColCount1);

			m_ScoreMx[MSAIndex1][MSAIndex2] = Score;
			m_ScoreMx[MSAIndex2][MSAIndex1] = Score;

			m_PathMx[MSAIndex1][MSAIndex2] = Path;
			m_PathMx[MSAIndex2][MSAIndex1] = InvertedPath;
			}
		}
	}

void PProg::LogPending(const string &s) const
	{
	Log("\nLogPending(%s) m_JoinIndex=%u\n", s.c_str(), m_JoinIndex);
	for (uint i = 0; i < SIZE(m_Pending); ++i)
		{
		uint Index = m_Pending[i];
		const MultiSequence &MSA = GetMSA(Index);
		uint SeqCount = MSA.GetSeqCount();
		uint ColCount = MSA.GetColCount();
		Log("  [%4u]  seqs=%u,cols=%u %s\n",
		  Index, SeqCount, ColCount, GetMSALabel(Index).c_str());
		}
	}

void PProg::Join_ByPrecomputedPath(uint Index1, uint Index2)
	{
	LogPending("Join start");

	asserta(SIZE(m_JoinMSAIndexes1) == m_JoinIndex);
	asserta(SIZE(m_JoinMSAIndexes2) == m_JoinIndex);
	m_JoinMSAIndexes1.push_back(Index1);
	m_JoinMSAIndexes2.push_back(Index2);

	uint NewMSAIndex = m_InputMSACount + m_JoinIndex;

	string NewMSALabel;
	Ps(NewMSALabel, "Join%u", m_JoinIndex+1);

	const string &MSALabel1 = m_MSALabels[Index1];
	const string &MSALabel2 = m_MSALabels[Index2];

	const MultiSequence &MSA1 = GetMSA(Index1);
	const MultiSequence &MSA2 = GetMSA(Index2);
	MultiSequence *MSA12 = new MultiSequence;
	const string &Path = m_PathMx[Index1][Index2];
	AlignMSAsByPath(MSA1, MSA2, Path, *MSA12);
	//AssertSeqsEq(MSA1, *MSA12);
	//AssertSeqsEq(MSA2, *MSA12);

	ProgressLog("Join %u/%u best pair %u, %u\n",
	  m_JoinIndex+1, m_JoinCount, Index1, Index2);
	Log("  Join_%u.X=%s\n", m_JoinIndex+1, MSALabel1.c_str());
	Log("  Join_%u.Y=%s\n", m_JoinIndex+1, MSALabel2.c_str());

	SetMSA(NewMSAIndex, *MSA12);
	SetMSALabel(NewMSAIndex, NewMSALabel);

	string JoinFileName = ".";
	if (optset_savedir)
		{
		string Prefix = opt(savedir);
		Dirize(Prefix);
		Ps(JoinFileName, "%sjoin%u", Prefix.c_str(), m_JoinIndex);
		ProgressLog("Writing join MSA: %s\n", JoinFileName.c_str());
		MSA12->WriteMFA(JoinFileName);
		}

	const uint JoinedSeqCount = MSA12->GetSeqCount();
	const uint JoinedColCount = MSA12->GetColCount();

	m_Pending.push_back(NewMSAIndex);

	uint PendingCountBeforeJoin = SIZE(m_Pending);
	DeleteIndexesFromPending(Index1, Index2);
	uint PendingCountAfterJoin = SIZE(m_Pending);
	asserta(PendingCountAfterJoin + 2 == PendingCountBeforeJoin);
	LogPending("Join end");
	}

void PProg::AlignNewToPending()
	{
	const uint PendingCount = SIZE(m_Pending);
	LogPending("AlignNewToPending");

	asserta(PendingCount > 0);
	const uint NewIndex = m_Pending[PendingCount-1];

	const string &NewMSALabel = m_MSALabels[NewIndex];
	const MultiSequence *NewMSA = &GetMSA(NewIndex);
	for (uint i = 0; i + 1 < PendingCount; ++i)
		{
		ProgressLog("Join %u/%u new vs. pending %u/%u\n\n",
		  m_JoinIndex+1, m_JoinCount, i+1, PendingCount);
		uint Index = m_Pending[i];
		const MultiSequence *MSA = &GetMSA(Index);

		const string &MSALabeli = m_MSALabels[i];
		string Path;
		float Score = AlignMSAsFlat(NewMSALabel + "+" + MSALabeli,
		  *NewMSA, *MSA, m_TargetPairCount, Path);

		string InvertedPath;
		InvertPath(Path, InvertedPath);

		m_ScoreMx[NewIndex][Index] = Score;
		m_ScoreMx[Index][NewIndex] = Score;

		m_PathMx[NewIndex][Index] = Path;
		m_PathMx[Index][NewIndex] = InvertedPath;
		}
	}

void PProg::WriteGuideTree(const string &FileName) const
	{
	if (FileName.empty())
		return;
	Tree GuideTree;
	MakeGuideTreeFromJoinOrder(m_JoinMSAIndexes1, m_JoinMSAIndexes2,
	  m_MSALabelToIndex, GuideTree);
	GuideTree.ToFile(FileName);
	}

//void PProg::CalcFwdFlat_PProg(uint GSI1, uint L1, 
//  uint GSI2, uint L2, float *Flat)
//	{
//	const byte *ByteSeq1 = GetByteSeqByGSI(GSI1);
//	const byte *ByteSeq2 = GetByteSeqByGSI(GSI2);
//	CalcFwdFlat(ByteSeq1, L1, ByteSeq2, L2, Flat);
//	}
//
//void PProg::CalcBwdFlat_PProg(uint GSI1, uint L1, 
//  uint GSI2, uint L2, float *Flat)
//	{
//	const byte *ByteSeq1 = GetByteSeqByGSI(GSI1);
//	const byte *ByteSeq2 = GetByteSeqByGSI(GSI2);
//	CalcBwdFlat(ByteSeq1, L1, ByteSeq2, L2, Flat);
//	}

void cmd_pprog()
	{
	PProg PP;
	vector<string> MSAFileNames;
	ReadStringsFromFile(opt(pprog), MSAFileNames);

	const uint MSACount = SIZE(MSAFileNames);
	asserta(MSACount > 1);
	const string &OutputFileName = opt(output);

	PP.m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	if (optset_paircount)
		PP.m_TargetPairCount = int(opt(paircount));

	bool IsNucleo;
	PP.LoadMSAs(MSAFileNames, IsNucleo);
	SetAlpha(IsNucleo ? ALPHA_Nucleo : ALPHA_Amino);

	InitProbcons();

	PP.Run();

	asserta(SIZE(PP.m_Pending) == 1);
	uint Index = PP.m_Pending[0];
	const MultiSequence &FinalMSA = PP.GetFinalMSA();
	FinalMSA.WriteMFA(OutputFileName);
	PP.WriteGuideTree(opt(guidetreeout));
	}
