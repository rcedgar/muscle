#include "muscle.h"
#include "mpc.h"

/////////////////////////////////////////////////////////////////
// One round of the consistency transformation:
//                     1
//    P'(x[i]-y[j]) = ---  sum   sum P(x[i]-z[k]) P(z[k]-y[j])
//                    |S| z in S  k
//
// where S = {x, y, all other input seqs}
/////////////////////////////////////////////////////////////////
void MPC::DoRelax(uint Iter)
	{
	const uint SeqCount = GetSeqCount();
	vector<vector<SparseMatrix*> > NewSparseMatrices(SeqCount, vector<SparseMatrix*>(SeqCount, NULL));

	const uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairIndex = 0;
	omp_lock_t Lock;
	omp_init_lock(&Lock);
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		const pair<uint, uint>& Pair = m_Pairs[PairIndex];
		uint i = Pair.first;
		uint j = Pair.second;

		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount, "Consistency(%u/%u)",
			Iter+1, m_ConsistencyIterCount);
		omp_unset_lock(&Lock);

		Sequence* seq1 = m_InputSeqs->GetSequence(i);
		Sequence* seq2 = m_InputSeqs->GetSequence(j);

	// get the original posterior matrix
		vector<float>* posteriorPtr = m_SparseMatrices[i][j]->GetPosterior();
		asserta(posteriorPtr != 0);
		vector<float>& posterior = *posteriorPtr;

		const uint L1 = seq1->GetLength();
		const uint L2 = seq2->GetLength();

	// contribution from the summation where z = x and z = y
		for (uint k = 0; k < (L1 + 1) * (L2 + 1); k++)
			posterior[k] += posterior[k];

	// contribution from all other m_InputSeqs
		for (uint k = 0; k < SeqCount; k++)
			{
			if (k != i && k != j)
				{
				if (k < i)
					Relax1(m_SparseMatrices[k][i], m_SparseMatrices[k][j], posterior);
				else if (k > i && k < j)
					Relax(m_SparseMatrices[i][k], m_SparseMatrices[k][j], posterior);
				else
					{
					SparseMatrix *temp =
						m_SparseMatrices[j][k]->ComputeTranspose();
					Relax(m_SparseMatrices[i][k], temp, posterior);
					delete temp;
					}
				}
			}

	// Normalize
		for (uint k = 0; k < (L1 + 1) * (L2 + 1); k++)
			posterior[k] /= SeqCount;

	// Delete cells not in original matrix
		SparseMatrix* matXY = m_SparseMatrices[i][j];
		for (uint y = 0; y <= L2; y++)
			posterior[y] = 0;

		for (uint x = 1; x <= L1; x++)
			{
			vector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
			vector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
			vector<float>::iterator base = posterior.begin() + x * (L2 + 1);
			uint curr = 0;
			while (XYptr != XYend)
				{
			// zero out all cells until the first filled column
				while (curr < (uint) XYptr->first)
					{
					base[curr] = 0;
					curr++;
					}

			// Skip this column
				curr++;
				++XYptr;
				}

		// zero out cells after last column
			while (curr <= L2)
				{
				base[curr] = 0;
				curr++;
				}
			}

		NewSparseMatrices[i][j] = new SparseMatrix(L1, L2, posterior);
		NewSparseMatrices[j][i] = NULL;

		delete posteriorPtr;
		}

// Replace previous posterior matrices
	for (uint i = 0; i < SeqCount; ++i)
		{
		for (uint j = 0; j < SeqCount; ++j)
			{
			delete m_SparseMatrices[i][j];
			m_SparseMatrices[i][j] = NewSparseMatrices[i][j];
			}
		}
	}
