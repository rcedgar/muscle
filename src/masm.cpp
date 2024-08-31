#include "muscle.h"
#include "masm.h"
#include "mega.h"

static float ScorePP(const MASMCol &PPA, const vector<byte> &ProfColB)
	{
	const uint FeatureCount = Mega::GetFeatureCount();
	assert(SIZE(ProfColB) == FeatureCount);

	float TotalScore = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		const vector<float> &ScoresA = PPA.m_ScoresVec[FeatureIdx];
		byte LetterB = ProfColB[FeatureIdx];
		TotalScore += ScoresA[LetterB];
		}
	
	return TotalScore;
	}

void MASM::MakeSMx(const vector<vector<byte> > &ProfB, Mx<float> &SMx) const
	{
	const uint LA = m_ColCount;
	const uint LB = SIZE(ProfB);
	SMx.AllocData(LA, LB);
	for (uint PosA = 0; PosA < LA; ++PosA)
		{
		const MASMCol &ColA = GetCol(PosA);
		for (uint PosB = 0; PosB < LB; ++PosB)
			{
			const vector<byte> &ColB = ProfB[PosB];
			float Score = ScorePP(ColA, ColB);
			SMx.Put(PosA, PosB, Score);
			}
		}
	}

void MASM::GetCounts(uint ColIndex, uint &LetterCount,
  uint &GapOpenCount, uint &GapExtCount, uint &GapCloseCount)
	{
	LetterCount = 0;
	GapOpenCount = 0;
	GapExtCount = 0;
	GapCloseCount = 0;

	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const char *Seq = m_Aln->GetCharPtr(SeqIndex);

		bool gap_prev = (ColIndex == 0 ? false : isgap(Seq[ColIndex-1]));
		bool gap = isgap(Seq[ColIndex]);
		bool gap_next = (ColIndex + 1 == m_ColCount ? false : isgap(Seq[ColIndex+1]));

		if (gap)
			{
			if (gap_prev)
				++GapExtCount;
			else if (gap_next)
				++GapOpenCount;
			else
				++GapCloseCount;
			}
		else
			++LetterCount;
		}
	}

void MASM::GetFreqs(uint ColIndex, uint FeatureIdx, vector<float> &Freqs)
	{
	Freqs.clear();
	const uint AlphaSize = Mega::GetAlphaSize(FeatureIdx);
	vector<uint> Counts(AlphaSize);

	for (uint SeqIdx = 0; SeqIdx < m_SeqCount; ++SeqIdx)
		{
		byte Letter = m_FeatureAlnVec[FeatureIdx][SeqIdx][ColIndex];
		if (Letter == UINT8_MAX)
			continue;
		asserta(Letter < AlphaSize);
		Counts[Letter] += 1;
		}

	for (byte Letter = 0; Letter < AlphaSize; ++Letter)
		{
		uint n = Counts[Letter];
		float Freq = float(n)/m_SeqCount;
		Freqs.push_back(Freq);
		}
	}

void MASM::GetFreqsVec(uint ColIndex, vector<vector<float> > &FreqsVec)
	{
	const uint FeatureCount = Mega::GetFeatureCount();
	FreqsVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		GetFreqs(ColIndex, FeatureIdx, FreqsVec[FeatureIdx]);
	}

void MASM::FromMSA(const MultiSequence &Aln, float GapOpen, float GapExt)
	{
	asserta(GapOpen >= 0);
	asserta(GapExt >= 0);

	Clear();

	m_Aln = &Aln;
	m_ColCount = Aln.GetColCount();
	m_SeqCount = Aln.GetSeqCount();
	m_FeatureCount = Mega::GetFeatureCount();
	m_AAFeatureIdx = Mega::GetAAFeatureIdx();
	SetUngappedSeqs();
	SetFeatureAlnVec();
	for (uint ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		uint LetterCount, GapOpenCount, GapExtCount, GapCloseCount;
		GetCounts(ColIndex, LetterCount, GapOpenCount, GapExtCount, GapCloseCount);
		asserta(LetterCount + GapOpenCount + GapExtCount + GapCloseCount == m_SeqCount);

		MASMCol &Col = *new MASMCol;
		Col.m_MASM = this;

		Col.m_LetterFreq = float(LetterCount)/m_SeqCount;
		Col.m_GapOpenFreq = float(GapOpenCount)/m_SeqCount;
		Col.m_GapExtFreq = float(GapExtCount)/m_SeqCount;
		Col.m_GapCloseFreq = float(GapCloseCount)/m_SeqCount;
		asserta(feq(Col.m_GapOpenFreq +
				    Col.m_GapExtFreq +
					Col.m_GapCloseFreq +
					Col.m_LetterFreq, 1));

		Col.m_GapOpen = (1 - Col.m_GapOpenFreq)*GapOpen/2;
		Col.m_GapClose = (1 - Col.m_GapCloseFreq)*GapOpen/2;
		Col.m_GapExt = (1 - Col.m_GapExtFreq)*GapExt/2;

		GetFreqsVec(ColIndex, Col.m_FreqsVec);
		Col.SetScoreVec();

		m_Cols.push_back(&Col);
		}
	}

void MASM::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToFile(f);
	CloseStdioFile(f);
	}

void MASM::ToFile(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "MSAM\t%u\t%u\t%u\n",
	  m_SeqCount, m_ColCount, Mega::GetFeatureCount());
	for (uint i = 0; i < m_ColCount; ++i)
		m_Cols[i]->ToFile(f, i);
	}

void MASM::SetFeatureAlnVec()
	{
	const uint FeatureCount = Mega::GetFeatureCount();
	m_FeatureAlnVec.clear();
	m_FeatureAlnVec.resize(FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		SetFeatureAln(FeatureIdx);
	}

void MASM::SetUngappedSeqs()
	{
	m_UngappedSeqs.clear();
	for (uint SeqIdx = 0; SeqIdx < m_SeqCount; ++SeqIdx)
		{
		string UngappedSeq;
		const Sequence &seq = *m_Aln->GetSequence(SeqIdx);
		asserta(seq.GetLength() == m_ColCount);
		const char *charseq = seq.GetCharPtr();
		for (uint Col = 0; Col < m_ColCount; ++Col)
			{
			char c = charseq[Col];
			if (!isgap(c))
				UngappedSeq += c;
			}
		m_UngappedSeqs.push_back(UngappedSeq);
		}
	}

void MASM::SetFeatureAln(uint FeatureIdx)
	{
	asserta(FeatureIdx < SIZE(m_FeatureAlnVec));
	vector<vector<byte> > &FeatureAln = m_FeatureAlnVec[FeatureIdx];
	FeatureAln.resize(m_SeqCount);
	asserta(SIZE(m_UngappedSeqs) == m_SeqCount);
	for (uint SeqIdx = 0; SeqIdx < m_SeqCount; ++SeqIdx)
		{
		const string &UngappedSeq = m_UngappedSeqs[SeqIdx];
		const vector<vector<byte> > *ptrMegaProfile = 
		  Mega::GetProfileBySeq(UngappedSeq, true);
		const vector<vector<byte> > &MegaProfile = *ptrMegaProfile;
		vector<byte> &Row = FeatureAln[SeqIdx];
		const Sequence &seq = *m_Aln->GetSequence(SeqIdx);
		asserta(seq.GetLength() == m_ColCount);
		const char *charseq = seq.GetCharPtr();
		uint Pos = 0;
		for (uint Col = 0; Col < m_ColCount; ++Col)
			{
			char c = charseq[Col];
			if (isgap(c))
				Row.push_back(UINT8_MAX);
			else
				{
				byte Letter = MegaProfile[Pos][FeatureIdx];
				Row.push_back(Letter);
				++Pos;
				}
			}
		}
	}

void MASM::MakeSMx_Sequence(const Sequence &Q, Mx<float> &SMx) const
	{
	const uint LM = m_ColCount;
	const uint LQ = Q.GetLength();
	SMx.Alloc(LM, LQ);
	const char *qs = Q.GetCharPtr();
	const uint FeatureCount = Mega::GetFeatureCount();
	for (uint Col = 0; Col < LM; ++Col)
		{
		const MASMCol &MC = *m_Cols[Col];
		const vector<float> &Scores = MC.GetAAScores();
		for (uint PosQ = 0; PosQ < LQ; ++PosQ)
			{
			char q = qs[PosQ];
			byte Letter = g_CharToLetterAmino[q];
			float Score = 0;
			if (Letter < 20)
				Score = Scores[Letter];
			SMx.Put(Col, PosQ, Score);
			}
		}
	}