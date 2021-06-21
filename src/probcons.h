#pragma once

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

#include "multisequence.h"
#include "ScoreType.h"
#include "pairhmm.h"
#include "guidetree.h"
#include "sparsematrix.h"

MultiSequence* RunMPC(MultiSequence *InputSeqs);
MultiSequence* AlignAlignments(MultiSequence* align1, MultiSequence* align2,
	const vector<vector<SparseMatrix*> >& sparseMatrices);
void Relax(SparseMatrix* matXZ, SparseMatrix* matZY, vector<float>& posterior);
void Relax1(SparseMatrix* matXZ, SparseMatrix* matZY, vector<float>& posterior);

void DoIterativeRefinement(const vector<vector<SparseMatrix*> >& sparseMatrices,
	MultiSequence*& alignment);
void InitProbcons();

void Relax(SparseMatrix* matXZ, SparseMatrix* matZY, vector<float>& posterior);
void Relax1(SparseMatrix* matZX, SparseMatrix* matZY, vector<float>& posterior);
void DoIterativeRefinement(const vector<vector<SparseMatrix*> >& sparseMatrices,
	MultiSequence*& alignment);

float AlignMSAs(const string &ProgressStr,
  const MultiSequence &MSA1, 
  const MultiSequence &MSA2,
  uint TargetPairCount,
  vector<char> &Path);

void ValidatePath(const vector<char> &Path, uint LX, uint LY);
void InvertPath(const vector<char> &Path, vector<char> &InvertedPath);
void InsertGappyPositions(const vector<char> &OccPath,
  uint FullColCount1, uint FullColCount2,
  vector<uint> &OccCols1, vector<uint> &OccCols2,
  vector<char> &FullPath);
void AlignMSAsByPath(const MultiSequence &MSA1, const MultiSequence &MSA2,
  const vector<char> &Path, MultiSequence &MSA12);
void _AssertSeqsEq(const char *FileName, uint LineNr,
  const MultiSequence &MSA1, const MultiSequence &MSA2);
#define AssertSeqsEq(MSA1, MSA2)	_AssertSeqsEq(__FILE__, __LINE__, MSA1, MSA2)

static const uint DEFAULT_CONSISTENCY_ITERS = 2;
static const uint DEFAULT_REFINE_ITERS = 100;
