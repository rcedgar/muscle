#pragma once

const float MIN_SPARSE_VALUE = 0.01f;

class SparseMx
    {
public:
    uint m_RowCount = 0;
    uint m_ColCount = 0;

    vector<int> m_RowSizes;
    vector<int> m_RowOffsets;
    vector<float> m_Values;
    vector<uint> m_ColIndexes;

public:
    SparseMx()
        {
        m_RowCount = 0;
        m_ColCount = 0;
        }

    ~SparseMx()
        {
        Clear();
        }

    void Clear()
        {
        m_RowCount = 0;
        m_ColCount = 0;

        m_RowSizes.clear();
        m_RowOffsets.clear();
        m_Values.clear();
        m_ColIndexes.clear();
        }

    void FromTomPosterior(const vector<float> &TomPosterior,
      uint L1, uint L2);

public:
    float GetValue(uint RowIndex, uint ColIndex) const;
    void ComputeTranspose(SparseMx &T) const;
    void GetTomPosterior(vector<float> &TomPosterior) const;
    };
