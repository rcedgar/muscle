#include "muscle.h"
#include "pathscorer.h"

float PathScorer::GetGlobalScore(uint LA, uint LB, const string &Path)
	{
	const uint ColCount = SIZE(Path);
	asserta(ColCount > 0);
	float Total = 0;
	uint PosA = 0;
	uint PosB = 0;
	char LastState = 'S';
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char State = Path[Col];
		Total += GetScore(LastState, State, PosA, PosB);
		switch (State)
			{
		case 'M':	++PosA; ++PosB; break;
		case 'D':	++PosA; break;
		case 'I':	++PosB; break;
		default:	asserta(false);
			}
		LastState = State;
		}
	asserta(PosA == LA && PosB == LB);
	Total += GetScore(LastState, 'E', PosA, PosB);
	return Total;
	}

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
		Total += GetScore(LastState, State, PosA, PosB);
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

void PathScorer::TermizePath(string &Path)
	{
	const uint ColCount = SIZE(Path);
	uint FirstM = UINT_MAX;
	uint LastM = UINT_MAX;
	for (uint i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = i;
			LastM = i;
			}
		}
	if (FirstM == UINT_MAX)
		return;
	for (uint Col = 0; Col < FirstM; ++Col)
		{
		char c = Path[Col];
		if (c == 'D')
			Path[Col] = 'd';
		else if (c == 'I')
			Path[Col] = 'i';
		}
	for (uint Col = LastM+1; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'D')
			Path[Col] = 'd';
		else if (c == 'I')
			Path[Col] = 'i';
		}
	}

float PathScorer::GetScore(char FromState, char ToState,
  uint PosA, uint PosB)
	{
	if (0)
		;

#define c(x, y)	else if (FromState == #x[0] && ToState == #y[0]) return GetScore##x##y(PosA, PosB);
	c(S, M)
	c(S, d)
	c(S, i)

	c(M, M)
	c(M, D)
	c(M, d)
	c(M, I)
	c(M, i)

	c(D, M)
	c(d, M)
	c(D, D)
	c(d, d)

	c(I, M)
	c(i, M)
	c(I, I)
	c(i, i)

	c(M, E)
	c(d, E)
	c(i, E)
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
	}

float PathScorer_MASM_Mega::GetScoreSM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_MASM_Mega::GetScoreSd(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreSi(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreMM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	return m;
	}

float PathScorer_MASM_Mega::GetScoreMD(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreMd(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreMI(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreMi(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreDM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	}

float PathScorer_MASM_Mega::GetScoredM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	}

float PathScorer_MASM_Mega::GetScoreDD(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoredd(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreIM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	}

float PathScorer_MASM_Mega::GetScoreiM(uint PosA, uint PosB)
	{
	float m = GetMatchScore(PosA, PosB);
	}

float PathScorer_MASM_Mega::GetScoreII(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreii(uint PosA, uint PosB)
	{
	}

float PathScorer_MASM_Mega::GetScoreME(uint PosA, uint PosB)
	{
	return 0;
	}

float PathScorer_MASM_Mega::GetScoredE(uint PosA, uint PosB)
	{
	return 0;
	}

float PathScorer_MASM_Mega::GetScoreiE(uint PosA, uint PosB)
	{
	return 0;
	}
