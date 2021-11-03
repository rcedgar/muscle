#pragma once

const float MIN_SPARSE_PROB = 0.01f;
const float MIN_SPARSE_SCORE = logf(MIN_SPARSE_PROB); // -4.6

class MySparseMx
	{
public:
	uint m_LX = 0;
	uint m_LY = 0;
	uint m_VecSize = 0;

	uint m_MaxVecSize = 0;
	byte *m_ValueVec = 0;

	uint m_MaxLX = 0;
	uint *m_Offsets = 0;

	const byte *m_X = 0;
	const byte *m_Y = 0;

public:
	MySparseMx()
		{
		m_LX = 0;
		m_LY = 0;
		m_VecSize = 0;

		m_MaxVecSize = 0;
		m_ValueVec = 0;

		m_MaxLX = 0;
		m_Offsets = 0;

		m_X = 0;
		m_Y = 0;
		}
	~MySparseMx()
		{
		Clear();
		}

	void Clear()
		{
		myfree(m_ValueVec);
		myfree(m_Offsets);
		m_ValueVec = 0;
		m_Offsets = 0;
		m_MaxVecSize = 0;
		m_MaxLX = 0;
		m_VecSize = 0;
		m_LX = 0;
		m_LY = 0;
		}

	float GetProb_Offset(uint Offset) const
		{
		const float *ptr_P = (float *) (m_ValueVec + 8*Offset);
		float P = *ptr_P;
		return P;
		}

	uint GetCol_Offset(uint Offset) const
		{
		const uint *ptr_j = (uint *) (m_ValueVec + 8*Offset + 4);
		uint Col = *ptr_j;
		return Col;
		}

	void SetProb_Offset(uint Offset, float P)
		{
		float *ptr_P = (float *) (m_ValueVec + 8*Offset);
		*ptr_P = P;
		}

	void SetCol_Offset(uint Offset, uint Col) const
		{
		uint *ptr_j = (uint *) (m_ValueVec + 8*Offset + 4);
		*ptr_j = Col;
		}

	void AllocLX(uint LX);
	void AllocVec(uint Size);
	void FromPost(const float *Post, uint LX, uint LY);
	void UpdateFromPost(const MySparseMx &OldMx,
	  const float *Post, uint SeqCount);
	void GetColToRowLoHi(vector<uint> &ColToRowLo, vector<uint> &ColToRowHi) const;
	void ToPost(float *Post) const;
	float GetProb(uint i, uint j) const;
	uint GetOffset(uint i) const;
	uint GetSize(uint i) const;
	float GetMaxProbRow(uint i) const;

    void LogMe() const;
	void LogStats(const char *Msg = "") const;
	const uint GetLX() const { return m_LX; }
	const uint GetLY() const { return m_LY; }
	};
