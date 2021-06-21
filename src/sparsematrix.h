#pragma once

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

const float POSTERIOR_CUTOFF = 0.01f;

// Sparse matrix entry, first=column, second=value
typedef pair<int,float> PIF;

class SparseMatrix
    {
    int seq1Length = 0;
    int seq2Length = 0;
    vector<int> rowSize;
    vector<PIF> data;
    vector<vector<PIF>::iterator> rowPtrs;

private:
  SparseMatrix (){}

public:
    SparseMatrix (int a_seq1Length, int a_seq2Length,
      const vector<float> &posterior);

    vector<PIF>::iterator GetRowPtr (int row) const;
    float GetValue (int row, int col) const;
    int GetRowSize (int row) const;
    int GetSeq1Length () const {return seq1Length;}
    int GetSeq2Length () const {return seq2Length;}
    int GetNumCells () const {return int(data.size());}
    void Print (ostream &outfile) const;
    SparseMatrix *ComputeTranspose () const;
    vector<float> *GetPosterior () const;

    void LogMe() const;
    };

static inline int PosteriorCoordsToOffset(int Pos1, int L1, int Pos2, int L2)
	{
	assert(Pos1 > 0 && Pos1 <= L1);
	assert(Pos2 > 0 && Pos2 <= L2);
	int Offset = Pos1*(L2 + 1) + Pos2;
	return Offset;
	}

void LogPosterior(const vector<float> &Posterior, int L1, int L2);
