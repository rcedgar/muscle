#include "muscle.h"
#include "sparsemx.h"

// Tom's posterior format is triangular with 1-based coordinates
// and all zeros (ignored) in first row and first column.
void SparseMx::FromTomPosterior(const vector<float> &TomPosterior,
  uint L1, uint L2)
	{
	uint L1_1 = L1 + 1;
	uint L2_1 = L2 + 1;
	uint TomSize = L1_1*L2_1;
	asserta(SIZE(TomPosterior) == TomSize);

	m_RowCount = L1;
	m_ColCount = L2;
	m_RowSizes.resize(m_RowCount, UINT_MAX);
	m_RowOffsets.resize(m_RowCount, UINT_MAX);

// Count cells above threshold
// Pre-processing step to avoid memory
// thrashing by push_back()'s.
	uint SparseCount = 0;
	for (uint RowIndex = 1; RowIndex <= L1; ++RowIndex)
		{
		uint PosteriorOffset = RowIndex*L2_1;
		for (uint ColIndex = 1; ColIndex <= L2; ++ColIndex)
			{
			float Value = TomPosterior[PosteriorOffset++];
			if (Value >= MIN_SPARSE_VALUE)
				++SparseCount;
			}
		}
	m_Values.resize(SparseCount, 0);
	m_ColIndexes.resize(SparseCount, UINT_MAX);

	uint SparseOffset = 0;
	for (uint RowIndex = 1; RowIndex <= L1; ++RowIndex)
		{
		uint PosteriorOffset = RowIndex*L2_1;
		uint RowSize = 0;
		m_RowOffsets[RowIndex] = SparseOffset;
		for (uint ColIndex = 1; ColIndex <= L2; ++ColIndex)
			{
			float Value = TomPosterior[PosteriorOffset++];
			if (Value >= MIN_SPARSE_VALUE)
				{
				m_Values[SparseOffset] = Value;
				m_ColIndexes[SparseOffset] = ColIndex;
				++RowSize;
				++SparseOffset;
				}
			}
		m_RowSizes[RowIndex] = RowSize;
		}
	asserta(SparseOffset == SparseCount);
	}

float SparseMx::GetValue(uint RowIndex, uint ColIndex) const
	{
	assert(RowIndex < m_RowCount);
	assert(ColIndex < m_ColCount);
	uint Offset = m_RowOffsets[RowIndex];
	uint RowSize = m_RowSizes[RowIndex];
	const uint SparseCount = SIZE(m_Values);
	asserta(Offset + RowSize <= SparseCount);
	for (uint i = 0; i < RowSize; ++i)
		{
		uint Col = m_ColIndexes[Offset++];
		if (Col >= ColIndex)
			{
			if (Col == ColIndex)
				return m_Values[Offset];
			return 0;
			}
		}
	return 0;
	}
