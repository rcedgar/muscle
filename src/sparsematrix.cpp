#include "myutils.h"
#include "probcons.h"

/***
Derived from PROBCONS code by Chuong B. Do.
http://probcons.stanford.edu/download.html
doi: 10.1101/gr.2821705
***/

// Constructor.  Builds a sparse matrix from a posterior matrix.
// Note that the expected format for the posterior matrix is as
// a (seq1Length+1) x (seq2Length+1) matrix where the 0th row
// and 0th column are ignored (they should contain all zeroes).
// Data is stored in vector "data", "rowPtrs" is a vector
// with one iterator per row pointing to the start of the row,
// "rowSize" is a vector with row sizes.
SparseMatrix::SparseMatrix(int a_seq1Length, int a_seq2Length,
  const vector<float>& posterior)
	{
	seq1Length = a_seq1Length;
	seq2Length = a_seq2Length;
	int numCells = 0;

	assert(seq1Length > 0);
	assert(seq2Length > 0);

	// calculate memory required; count the number of cells in the
	// posterior matrix above the threshold
	vector<float>::const_iterator postPtr = posterior.begin();
	for (int i = 0; i <= seq1Length; i++)
		{
		for (int j = 0; j <= seq2Length; j++)
			{
			if (*(postPtr++) >= POSTERIOR_CUTOFF)
				{
				assert(i != 0 && j != 0);
				numCells++;
				}
			}
		}

	// allocate memory
	data.resize(numCells);

	rowSize.resize(seq1Length + 1);
	rowSize[0] = -1;

	rowPtrs.resize(seq1Length + 1); 
	rowPtrs[0] = data.end();

	// build sparse matrix
	// note that we're skipping the first row here...
	postPtr = posterior.begin() + seq2Length + 1;
	vector<PIF>::iterator dataPtr = data.begin();
	for (int i = 1; i <= seq1Length; i++)
		{
		// ...and skipping the first column of each row
		postPtr++;
		rowPtrs[i] = dataPtr;
		for (int j = 1; j <= seq2Length; j++)
			{
			if (*postPtr >= POSTERIOR_CUTOFF)
				{
				dataPtr->first = j;
				dataPtr->second = *postPtr;
				dataPtr++;
				}
			postPtr++;
			}
		rowSize[i] = int(dataPtr - rowPtrs[i]);
		}
	}

/////////////////////////////////////////////////////////////////
// SparseMatrix::GetValue()
//
// Returns value at a particular row, column.
/////////////////////////////////////////////////////////////////
float SparseMatrix::GetValue(int row, int col) const
	{
	assert(row >= 1 && row <= seq1Length);
	assert(col >= 1 && col <= seq2Length);
	for (int i = 0; i < rowSize[row]; i++)
		{
		if (rowPtrs[row][i].first == col)
			return rowPtrs[row][i].second;
		}
	return 0;
	}

/////////////////////////////////////////////////////////////////
// SparseMatrix::GetRowSize()
//
// Returns the number of entries in a particular row.
/////////////////////////////////////////////////////////////////

int SparseMatrix::GetRowSize(int row) const
	{
	assert(row >= 1 && row <= seq1Length);
	return rowSize[row];
	}


void SparseMatrix::Print(ostream& outfile) const
	{
	outfile << "Sparse Matrix:" << endl;
	for (int i = 1; i <= seq1Length; i++)
		{
		outfile << "  " << i << ":";
		for (int j = 0; j < rowSize[i]; j++)
			{
			outfile << " (" << rowPtrs[i][j].first << ","
				<< rowPtrs[i][j].second << ")";
			}
		outfile << endl;
		}
	}

/////////////////////////////////////////////////////////////////
// SparseMatrix::ComputeTranspose()
//
// Returns a new sparse matrix containing the transpose of the
// current matrix.
/////////////////////////////////////////////////////////////////

SparseMatrix* SparseMatrix::ComputeTranspose() const
	{
	// create a new sparse matrix
	SparseMatrix* ret = new SparseMatrix();
	int numCells = int(data.size());

	ret->seq1Length = seq2Length;
	ret->seq2Length = seq1Length;

	// allocate memory
	ret->data.resize(numCells);
	ret->rowSize.resize(seq2Length + 1); ret->rowSize[0] = -1;
	ret->rowPtrs.resize(seq2Length + 1); ret->rowPtrs[0] = ret->data.end();

	// compute row sizes
	for (int i = 1; i <= seq2Length; i++) ret->rowSize[i] = 0;
	for (int i = 0; i < numCells; i++)
		ret->rowSize[data[i].first]++;

	// compute row ptrs
	for (int i = 1; i <= seq2Length; i++) {
		ret->rowPtrs[i] = (i == 1) ? ret->data.begin() : ret->rowPtrs[i - 1] + ret->rowSize[i - 1];
		}

	// now fill in data
	vector<vector<PIF>::iterator> currPtrs = ret->rowPtrs;

	for (int i = 1; i <= seq1Length; i++) {
		vector<PIF>::iterator row = rowPtrs[i];
		for (int j = 0; j < rowSize[i]; j++) {
			currPtrs[row[j].first]->first = i;
			currPtrs[row[j].first]->second = row[j].second;
			currPtrs[row[j].first]++;
			}
		}

	return ret;
	}

vector<float>* SparseMatrix::GetPosterior() const
	{
	// create a new posterior matrix
	vector<float>* posteriorPtr =
	  new vector<float>((seq1Length + 1) * (seq2Length + 1)); 
	assert(posteriorPtr);
	vector<float>& posterior = *posteriorPtr;

	// build the posterior matrix
	for (int i = 0; i < (seq1Length + 1) * (seq2Length + 1); i++)
		posterior[i] = 0;

	for (int i = 1; i <= seq1Length; i++)
		{
		vector<float>::iterator postPtr = 
		  posterior.begin() + i * (seq2Length + 1);

		for (int j = 0; j < rowSize[i]; j++)
			postPtr[rowPtrs[i][j].first] = rowPtrs[i][j].second;
		}
	return posteriorPtr;
	}

vector<PIF>::iterator SparseMatrix::GetRowPtr(int row) const
	{
	assert(row >= 1 && row <= seq1Length);
	return rowPtrs[row];
	}

void LogPosterior(const vector<float> &Posterior, int L1, int L2)
	{
	Log("\n");
	Log("Posterior(%p) L1 %u, L2 %u\n", &Posterior, L1, L2);
	Log("  Row");
	for (int i = 1; i <= L2; ++i)
		Log("  %8d", i);
	Log("\n");
	for (int Pos1 = 1; Pos1 <= L1; ++Pos1)
		{
		Log("%5d", Pos1);
		for (int Pos2 = 1; Pos2 <= L2; ++Pos2)
			{
			int Offset = PosteriorCoordsToOffset(Pos1, L1, Pos2, L2);
			float Value = Posterior[Offset];
			if (Value == 0)
				Log("  %8.8s", ".");
			else
				Log("  %8.3g", Value);
			}
		Log("\n");
		}
	}

void SparseMatrix::LogMe() const
	{
	Log("\n");
	Log("SparseMatrix(%p) L1=%u, L2=%u\n",
	  this, seq1Length, seq2Length);

	asserta(SIZE(rowSize) == seq1Length+1);
	asserta(SIZE(rowPtrs) == seq1Length+1);
	for (int RowIndex = 1; RowIndex <= seq1Length; ++RowIndex)
		{
		int Size = rowSize[RowIndex];
		vector<PIF>::const_iterator p = rowPtrs[RowIndex];

		Log("[Row %5d]", RowIndex);
		Log("  [size %5d]", Size);
		for (int i = 0; i < Size; ++i)
			{
			const PIF &pif = *p;
			Log(" %d=%.3g", pif.first, pif.second);
			}
		Log("\n");
		}
	Log("\n");
	Log("  Row   Size");
	for (int i = 1; i <= seq2Length; ++i)
		Log("  %8d", i);
	Log("\n");

	vector<float> &Pos = *GetPosterior();
	for (int Pos1 = 1; Pos1 <= seq1Length; ++Pos1)
		{
		int Size = rowSize[Pos1];
		vector<PIF>::const_iterator p = rowPtrs[Pos1];

		Log("%5d", Pos1);
		Log("  %5d", Size);
		for (int Pos2 = 1; Pos2 <= seq2Length; ++Pos2)
			{
			int Offset = PosteriorCoordsToOffset(Pos1, seq1Length, Pos2, seq2Length);
			float Value = Pos[Offset];
			if (Value == 0)
				Log("  %8.8s", ".");
			else
				Log("  %8.3g", Value);
			}
		Log("\n");
		}
	delete &Pos;
	}
