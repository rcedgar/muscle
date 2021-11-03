#include "myutils.h"
#include "mysparsemx.h"

void MySparseMx::AllocVec(uint Size)
	{
	if (Size <= m_MaxVecSize)
		return;
	if (m_MaxVecSize > 0)
		myfree(m_ValueVec);

	m_MaxVecSize = Size + 256;
	asserta(sizeof(float) == 4 && sizeof(uint) == 4);
	m_ValueVec = myalloc(byte, m_MaxVecSize*8);
	}

float MySparseMx::GetMaxProbRow(uint i) const
	{
	const uint Offset = GetOffset(i);
	const uint Size = GetSize(i);
	float Max = 0;
	for (uint k = 0; k < Size; ++k)
		{
		float P = GetProb_Offset(Offset + k);
		Max = max(Max, P);
		}
	return Max;
	}

uint MySparseMx::GetOffset(uint i) const
	{
	assert(i < m_LX);
	uint Offset = m_Offsets[i];
	return Offset;
	}

uint MySparseMx::GetSize(uint i) const
	{
	assert(i < m_LX);
	uint Offset = GetOffset(i);
	uint Size = m_Offsets[i+1] - Offset;
	return Size;
	}

float MySparseMx::GetProb(uint i, uint j) const
	{
	uint Offset = GetOffset(i);
	uint Size = GetSize(i);
	for (uint k = 0; k < Size; ++k)
		{
		const uint *ptr_j = (uint *) (m_ValueVec + 8*Offset + 4);
		uint j2 = *ptr_j;
		if (j2 == j)
			{
			const float *ptr_P = (float *) (m_ValueVec + 8*Offset);
			return *ptr_P;
			}
		else if (j2 > j)
			return 0;
		++Offset;
		}
	return 0;
	}

void MySparseMx::AllocLX(uint LX)
	{
#if 0//TRACE
	Log("%p->AllocLX(%u) max %u\n", this, LX, m_MaxLX);
#endif
	if (LX <= m_MaxLX)
		return;

	if (m_MaxLX > 0)
		myfree(m_Offsets);

	m_MaxLX = LX + 128;
	m_Offsets = myalloc(uint, m_MaxLX+1);
#if 0//TRACE
	Log("%p->AllocLX(%u) newmax %u m_Offsets=%p\n", this, LX, m_MaxLX, m_Offsets);
#endif

#if DEBUG
	for (uint i = 0; i <= m_MaxLX; ++i)
		m_Offsets[i] = UINT_MAX;
#endif
	}

void MySparseMx::UpdateFromPost(const MySparseMx &OldMx,
  const float *Post, uint SeqCount)
	{
	uint VecSize = OldMx.m_VecSize;
	uint LX = OldMx.GetLX();
	uint LY = OldMx.GetLY();
	AllocLX(LX);
	AllocVec(VecSize);
	m_LX = LX;
	m_LY = LY;
	for (uint i = 0; i < LX; ++i)
		m_Offsets[i] = OldMx.m_Offsets[i];
	m_Offsets[LX] = OldMx.m_Offsets[LX];

	for (uint i = 0; i < m_LX; ++i)
		{
		uint Offset = GetOffset(i);
		uint Size = GetSize(i);
		for (uint k = 0; k < Size; ++k)
			{
			uint Col = OldMx.GetCol_Offset(Offset + k);
			float P = Post[i*LY + Col]/SeqCount;
			SetProb_Offset(Offset + k, P);
			SetCol_Offset(Offset + k, Col);
			}
		}
	}

void MySparseMx::FromPost(const float *Post, uint LX, uint LY)
	{
	asserta(sizeof(float) == sizeof(uint));
	m_LX = LX;
	m_LY = LY;

	AllocLX(LX);
	uint Offset = 0;
	for (uint i = 0; i < LX; ++i)
		{
		m_Offsets[i] = Offset;
		for (uint j = 0; j < LY; ++j)
			if (Post[i*LY + j] >= MIN_SPARSE_PROB)
				++Offset;
		}
	m_Offsets[LX] = Offset;

	m_VecSize = Offset;
	AllocVec(m_VecSize);

	Offset = 0;
	for (uint i = 0; i < LX; ++i)
		{
		for (uint j = 0; j < LY; ++j)
			{
			float P = Post[i*LY + j];
			if (P >= MIN_SPARSE_PROB)
				{
				float *ptr_P = (float *) (m_ValueVec + 8*Offset);
				uint *ptr_j = (uint *) (m_ValueVec + 8*Offset + 4);
				*ptr_P = P;
				*ptr_j = j;
				++Offset;
				}
			}
		}
	asserta(Offset == m_VecSize);
	}

void MySparseMx::LogStats(const char *Msg) const
	{
	Log("MySparseMx(%s) LX=%u, LY=%u VecSize=%u\n",
	  Msg, m_LX, m_LY, m_VecSize);
	}

void MySparseMx::LogMe() const
	{
	Log("\n");
	Log("MySparseMx(%p) LX=%u, LY=%u\n", this, m_LX, m_LY);

#if 0
	for (uint i = 0; i < m_LX; ++i)
		{
		uint Offset = GetOffset(i);
		uint Size = GetSize(i);

		Log("[Row %5u]", i);
		Log("  [off %5u]", Offset);
		Log("  [size %5u]", Size);
		for (uint k = 0; k < Size; ++k)
			{
			const float *ptr_P = (float *) (m_ValueVec + 8*Offset);
			const uint *ptr_j = (uint *) (m_ValueVec + 8*Offset + 4);
			Log(" %u=%.3g", *ptr_j, *ptr_P);
			++Offset;
			}
		Log("\n");
		}
#endif

	Log("\n");
	Log("  Row   Size");
	if (m_X != 0)
		Log("  x  ");
	for (uint j = 0; j < m_LY; ++j)
		Log("  %8u", j);
	Log("\n");
	if (m_Y != 0)
		{
		Log("\n");
		Log("                 ");
		for (uint j = 0; j < m_LY; ++j)
			Log("  %8c", m_Y[j]);
		Log("\n");
		}

	for (uint i = 0; i < m_LX; ++i)
		{
		uint Size = GetSize(i);
		Log("%5u", i);
		Log("  %5u", Size);
		if (m_X != 0)
			Log("  %c  ", m_X[i]);
		for (uint j = 0; j < m_LY; ++j)
			{
			float P = GetProb(i, j);
			if (P == 0)
				Log("  %8.8s", ".");
			else
				Log("  %8.3g", P);
			}
		Log("\n");
		}
	}

void MySparseMx::ToPost(float *Post) const
	{
	for (uint i = 0; i < m_LX*m_LY; ++i)
		Post[i] = 0;

	for (uint Row = 0; Row < m_LX; ++Row)
		{
		uint Offset = GetOffset(Row);
		uint Size = GetSize(Row);
		for (uint k = 0; k < Size; ++k)
			{
			float P = GetProb_Offset(Offset + k);
			uint Col = GetCol_Offset(Offset + k);
			Post[Row*m_LY + Col] = P;
			}
		}
	}

void MySparseMx::GetColToRowLoHi(vector<uint> &ColToRowLo, 
  vector<uint> &ColToRowHi) const
	{
	ColToRowLo.clear();
	ColToRowHi.clear();
	ColToRowLo.resize(m_LY, UINT_MAX);
	ColToRowHi.resize(m_LY, UINT_MAX);

	for (uint Row = 0; Row < m_LX; ++Row)
		{
		uint Offset = GetOffset(Row);
		uint Size = GetSize(Row);
		for (uint k = 0; k < Size; ++k)
			{
			uint Col = GetCol_Offset(Offset + k);
			assert(Col < m_LY);
			if (ColToRowLo[Col] == UINT_MAX)
				{
				ColToRowLo[Col] = Row;
				ColToRowHi[Col] = Row;
				}
			else
				{
				if (Row < ColToRowLo[Col])
					ColToRowLo[Col] = Row;
				if (Row > ColToRowHi[Col])
					ColToRowHi[Col] = Row;
				}
			}
		}
	}
