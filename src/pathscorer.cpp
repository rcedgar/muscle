#include "muscle.h"
#include "pathscorer.h"

float PathScorer::GetLocalScore(uint PosA, uint PosB, uint LA, uint LB,
  const string &Path)
	{
	const uint ColCount = SIZE(Path);
	asserta(ColCount > 0);
	asserta(Path[0] == 'M');
	asserta(Path[ColCount-1] == 'M');

	float Total = 0;
	char LastState = 'M';
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char State = Path[Col];
		float Score = GetScore(LastState, State, PosA, PosB);
		Total += Score;
		if (m_Trace)
			Log("Col %3u  PosA %3u  PosB %3u  %c%c  %10.3g  %10.3g\n",
			  Col, PosA, PosB, LastState, State, Score, Total);
		switch (State)
			{
		case 'M':	++PosA; ++PosB; break;
		case 'D':	++PosA; break;
		case 'I':	++PosB; break;
		default:	asserta(false);
			}
		LastState = State;
		}
	asserta(PosA <= LA && PosB <= LB);
	return Total;
	}

float PathScorer::GetScore(char FromState, char ToState,
  uint PosA, uint PosB)
	{
	if (0)
		;

#define c(x, y)	else if (FromState == #x[0] && ToState == #y[0]) return GetScore##x##y(PosA, PosB);
	c(M, M)
	c(M, D)
	c(M, I)

	c(D, M)
	c(D, D)

	c(I, M)
	c(I, I)
#undef c

	else
		asserta(false);
	return FLT_MAX;
	}

// Count match score in GetScoreXM(), not GetScoreMX()
//    (either should work)
// D is gap in MegaProfile
// I is gap in MASM

float PathScorer_MASM_Mega::GetMatchScore(uint PosA, uint PosB)
	{
	asserta(PosA < m_MASM->GetColCount());
	asserta(PosB < SIZE(*m_MegaProfile));
	const MASMCol &MCol = m_MASM->GetCol(PosA);
	const vector<byte> &PPos = (*m_MegaProfile)[PosB];
	float Score = MCol.GetMatchScore_MegaProfilePos(PPos);
	return Score;
	}

float PathScorer_MASM_Mega::GetScoreMM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_MASM_Mega::GetScoreMD(uint PosA, uint PosB)
	{
	asserta(PosA < m_MASM->GetColCount());
	const MASMCol &MCol = m_MASM->GetCol(PosA);
	return MCol.m_GapOpen;
	}

float PathScorer_MASM_Mega::GetScoreMI(uint PosA, uint PosB)
	{
	float Score = m_MASM->m_GapOpen/2;
	asserta(Score != FLT_MAX);
	asserta(Score < 0);
	return Score;
	}

float PathScorer_MASM_Mega::GetScoreDM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	asserta(PosA < m_MASM->GetColCount());
	const MASMCol &MCol = m_MASM->GetCol(PosA);
	float Score = m + MCol.m_GapClose;
	return Score;
	}

float PathScorer_MASM_Mega::GetScoreDD(uint PosA, uint PosB)
	{
	asserta(PosA < m_MASM->GetColCount());
	const MASMCol &MCol = m_MASM->GetCol(PosA);
	return MCol.m_GapExt;
	}

float PathScorer_MASM_Mega::GetScoreIM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	float Score = m_MASM->m_GapOpen/2;
	asserta(Score != FLT_MAX);
	asserta(Score < 0);
	return m + Score;
	}

float PathScorer_MASM_Mega::GetScoreII(uint PosA, uint PosB)
	{
	float Score = m_MASM->m_GapExt;
	asserta(Score != FLT_MAX);
	asserta(Score < 0);
	return Score;
	}

float PathScorer_AA_BLOSUM62::GetMatchScore(uint PosA, uint PosB)
	{
	float GetBlosumScoreChars(byte a, byte b);
	asserta(PosA < SIZE(m_SeqA));
	asserta(PosB < SIZE(m_SeqB));
	byte a = (byte) m_SeqA[PosA];
	byte b = (byte) m_SeqB[PosB];
	float Score = GetBlosumScoreChars(a, b);
	if (m_Trace)
		Log("\n    B62[%c,%c]=%.3g\n", a, b, Score);
	return Score;
	}

float PathScorer_AA_BLOSUM62::GetScoreMM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_AA_BLOSUM62::GetScoreMD(uint PosA, uint PosB)
	{
	return m_GapOpen;
	}

float PathScorer_AA_BLOSUM62::GetScoreMI(uint PosA, uint PosB)
	{
	return m_GapOpen;
	}

float PathScorer_AA_BLOSUM62::GetScoreDM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_AA_BLOSUM62::GetScoreDD(uint PosA, uint PosB)
	{
	return m_GapExt;
	}

float PathScorer_AA_BLOSUM62::GetScoreIM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_AA_BLOSUM62::GetScoreII(uint PosA, uint PosB)
	{
	return m_GapExt;
	}
