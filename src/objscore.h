#ifndef ObjScore_h
#define ObjScore_h

SCORE ScoreSeqPairGaps(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2);
SCORE ScoreSeqPairLetters(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2);
SCORE ScoreGaps(const MSA &msa, const unsigned Cols[], unsigned ColCount);

SCORE ObjScore(const MSA &msa, const unsigned SeqIndexes1[],
  unsigned uSeqCount1, const unsigned SeqIndexes2[], unsigned uSeqCount2);

SCORE ObjScoreIds(const MSA &msa, const unsigned Ids1[],
  unsigned uCount1, const unsigned Ids2[], unsigned uCount2);

void GetLetterScores(const MSA &msa, SCORE LetterScores[]);

SCORE ObjScoreDP(const MSA &msa1, const MSA &msa2, SCORE MatchScore[] = 0);
SCORE ObjScorePS(const MSA &msa, SCORE MatchScore[] = 0);
SCORE ObjScoreSP(const MSA &msa, SCORE MatchScore[] = 0);
SCORE ObjScoreXP(const MSA &msa, const MSA &msa2);
SCORE ObjScoreSPDimer(const MSA &msa);
SCORE ObjScoreDP_Profs(const ProfPos *PA, const ProfPos *PB, unsigned uColCount,
  SCORE MatchScore[] = 0);

SCORE DiffObjScore(
  const MSA &msa1, const PWPath &Path1, const unsigned Edges1[], unsigned uEdgeCount1, 
  const MSA &msa2, const PWPath &Path2, const unsigned Edges2[], unsigned uEdgeCount2);

#endif // ObjScore_h
