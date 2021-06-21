#ifndef pairhmm_h
#define pairhmm_h

#include "scoretype.h"
#include "sparsematrix.h"
#include "multisequence.h"
#include "hmmparams.h"

using namespace std;

enum HMMSTATE
	{
	HMMSTATE_M = 0,
	HMMSTATE_ISX = 1,
	HMMSTATE_ISY = 2,
	HMMSTATE_ILX = 3,
	HMMSTATE_ILY = 4,
	HMMSTATE_COUNT = 5
	};

static const uint InsertStateCount = 2;

static const float WILDCARD_PROB = 1e-5f;

class PairHMM
	{
public:
	static float m_StartScore[HMMSTATE_COUNT];
	static float m_TransScore[HMMSTATE_COUNT][HMMSTATE_COUNT];
	static float m_MatchScore[256][256];
	static float m_InsScore[256];

public:
	static void Create(const vector<float>& initDistribMat,
	  const vector<float>& gapOpen, const vector<float>& gapExtend,
	  const vector<vector<float>>& emitPairs,
	  const vector<float>& emitSingle);

	static void Create2(const vector<float>& initDistribMat,
	  const vector<vector<float> > &transMat, const vector<float>& emitSingle,
	  const vector<vector<float> > &emitPairs);

	static vector<float>* ComputeForwardMatrix(Sequence* seq1, Sequence* seq2);

	static vector<float>* ComputeBackwardMatrix(Sequence* seq1, Sequence* seq2);

	static float ComputeTotalProbability(int seq1Length, int seq2Length,
		const vector<float>& forward, const vector<float>& backward);

	static vector<float>* ComputePosteriorMatrix(Sequence* seq1, Sequence* seq2,
		const vector<float>& forward, const vector<float>& backward);

	static pair<vector<char>*, float> ComputeAlignment(int seq1Length, int seq2Length,
		const vector<float>& posterior);

	static vector<float>* BuildPosterior(MultiSequence* align1, MultiSequence* align2,
		const vector<vector<SparseMatrix*> >& sparseMatrices);

	static vector<float>* BuildPosterior2(const MultiSequence* align1,
	  const MultiSequence* align2,
	  const vector<vector<SparseMatrix*> >& sparseMatrices);

	static vector<float>* BuildPosterior3(
	  const MultiSequence* align1,
	  const MultiSequence* align2,
	  const vector<int> &SeqIndexes1,
	  const vector<int> &SeqIndexes2,
	  const vector<SparseMatrix*>& sparseMatrices);

	static float AlignPair(Sequence *X, Sequence *Y, vector<char> *Path = 0);
	static float AlignPair_StrPath(Sequence *X, Sequence *Y, string &Path);

	static void WriteParamsReport(const string &FileName);
	static void WriteParamsReport(FILE *f);
	};

#endif
