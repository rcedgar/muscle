#include "muscle.h"
#include "mpc.h"

float MPC::CalcAccAlnRow(uint PairIndex)
	{
	asserta(PairIndex < SIZE(m_Pairs));
	asserta(PairIndex < SIZE(m_AccAlnRows));
	const pair<uint, uint> &Pair = m_Pairs[PairIndex];
	uint MSASeqIndex1 = Pair.first;
	uint MSASeqIndex2 = Pair.second;

	vector<float> &AccAlnRow = m_AccAlnRows[PairIndex];
	AccAlnRow.clear();

	const Sequence *AlignedSeq1 = m_MSA->GetSequence(MSASeqIndex1);
	const Sequence *AlignedSeq2 = m_MSA->GetSequence(MSASeqIndex2);
	
	uint SMI1 = AlignedSeq1->GetSMI();
	uint SMI2 = AlignedSeq2->GetSMI();

	const SparseMatrix *ptrMx = m_SparseMatrices[SMI1][SMI2];
	if (ptrMx == 0)
		{
		ptrMx = m_SparseMatrices[SMI2][SMI1];
		asserta(ptrMx != 0);

		swap(MSASeqIndex1, MSASeqIndex2);
		swap(SMI1, SMI2);
		swap(AlignedSeq1, AlignedSeq2);
		}

	const Sequence *UnalignedSeq1 = m_InputSeqs->GetSequence(SMI1);
	const Sequence *UnalignedSeq2 = m_InputSeqs->GetSequence(SMI2);

	const uint L1 = UnalignedSeq1->GetLength();
	const uint L2 = UnalignedSeq2->GetLength();

	const SparseMatrix &Mx = *ptrMx;
	asserta(Mx.GetSeq1Length() == L1);
	asserta(Mx.GetSeq2Length() == L2);

	asserta(m_MSA != 0);
	const uint ColCount = m_MSA->GetColCount();

	const byte *Row1 = AlignedSeq1->GetBytePtr();
	const byte *Row2 = AlignedSeq2->GetBytePtr();

// 1-based positions in sparse matrices
	uint Pos1 = 1;
	uint Pos2 = 1;
	float SumP = 0;
	uint N = 0;

	for (uint Col = 0; Col < ColCount; ++Col)
		{
		byte c1 = Row1[Col];
		byte c2 = Row2[Col];

		float P = 0;
		bool Gap1 = isgap(c1);
		bool Gap2 = isgap(c2);
		if (!Gap1 && !Gap2)
			{
			P = Mx.GetValue(Pos1, Pos2);
			SumP += P;
			++N;
			}
		AccAlnRow.push_back(P);

		if (!Gap1)
			++Pos1;
		if (!Gap2)
			++Pos2;
		}
	float AvgEA = 0;
	if (N > 0)
		AvgEA = SumP/N;
	return AvgEA;
	}

void MPC::CalcAccAlnRows()
	{
	m_EAs.clear();
	m_AccAlnRows.clear();
	const uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	m_AccAlnRows.resize(PairCount);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		float EA = CalcAccAlnRow(PairIndex);
		m_EAs.push_back(EA);
		}
	}

void MPC::WriteAccAln1(FILE *f, uint PairIndex) const
	{
	if (f == 0)
		return;

	const pair<uint, uint> &Pair = m_Pairs[PairIndex];
	uint MSASeqIndex1 = Pair.first;
	uint MSASeqIndex2 = Pair.second;

	const string &Label1 = m_MSA->GetLabel(MSASeqIndex1);
	const string &Label2 = m_MSA->GetLabel(MSASeqIndex2);
	const vector<float> &AccAlnRow = m_AccAlnRows[PairIndex];
	float EA = m_EAs[PairIndex];
	const Sequence *MSASeq1 = m_MSA->GetSequence(MSASeqIndex1);
	const Sequence *MSASeq2 = m_MSA->GetSequence(MSASeqIndex2);
	const uint ColCount = m_MSA->GetColCount();
	asserta(ColCount > 0);
	asserta(SIZE(AccAlnRow) == ColCount);
	uint SMI1 = MSASeq1->GetSMI();
	uint SMI2 = MSASeq2->GetSMI();
	const Sequence *InputSeq1 = m_InputSeqs->GetSequence(SMI1);
	const Sequence *InputSeq2 = m_InputSeqs->GetSequence(SMI2);
	asserta(InputSeq1->GetLabel() == Label1);
	asserta(InputSeq2->GetLabel() == Label2);
	uint L1 = InputSeq1->GetLength();
	uint L2 = InputSeq2->GetLength();
	const byte *Row1 = MSASeq1->GetBytePtr();
	const byte *Row2 = MSASeq2->GetBytePtr();

	string OutRow1;
	string OutRow2;
	string OutRowA1;
	string OutRowA2;
	uint IdCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c1 = toupper(Row1[Col]);
		char c2 = toupper(Row2[Col]);
		float EA = AccAlnRow[Col];
		char a1 = ' ';
		char a2 = ' ';
		if (isgap(c1) && isgap(c2))
			continue;

		if (!isgap(c1) && !isgap(c2))
			{
			if (c1 == c2)
				++IdCount;
			if (EA > 0.995)
				{
				a1 = '^';
				a2 = '^';
				}
			else
				{
				uint PctEA = uint(EA*100.0f + 0.5f);
				if (PctEA >= 100)
					PctEA = 99;
				uint d1 = PctEA/10;
				uint d2 = PctEA%10;
				asserta(d1 >= 0 && d1 <= 9);
				asserta(d2 >= 0 && d2 <= 9);
				a1 = "0123456789"[d1];
				a2 = "0123456789"[d2];
				}
			}
		OutRow1 += c1;
		OutRow2 += c2;
		OutRowA1 += a1;
		OutRowA2 += a2;
		}

	float PctId = IdCount*100.0f/ColCount;
	fprintf(f, "\n");
	fprintf(f, "%u	%s	%s	%.3g	%.1f%%\n",
	  PairIndex, Label1.c_str(), Label2.c_str(), EA, PctId);
	fprintf(f, "%u	%s\n", PairIndex, OutRow1.c_str());
	fprintf(f, "%u	%s\n", PairIndex, OutRow2.c_str());
	fprintf(f, "%u	%s\n", PairIndex, OutRowA1.c_str());
	fprintf(f, "%u	%s\n", PairIndex, OutRowA2.c_str());
	}

void MPC::WriteAccAln(const string &FileName) const
	{
	if (FileName.empty())
		return;

	FILE *f = CreateStdioFile(FileName);
	const uint PairCount = SIZE(m_EAs);
	asserta(SIZE(m_AccAlnRows) == PairCount);
	asserta(SIZE(m_EAs) == PairCount);

	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		WriteAccAln1(f, PairIndex);

	CloseStdioFile(f);
	}
