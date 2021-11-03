#include "myutils.h"
#include "alpha.h"
#include "msa.h"

void DeleteAllGapColumns(MSA &M)
	{
	const uint SeqCount = M.GetSeqCount();
	const uint ColCount = M.GetColCount();

	vector<bool> Keeps;
	uint KeepCount = 0;
	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		bool Keep = (M.GetGapCount(ColIndex) < SeqCount);
		Keeps.push_back(Keep);
		if (Keep)
			++KeepCount;
		}
	if (KeepCount == 0)
		Die("MSA is all gaps");

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char *Seq = M.m_szSeqs[SeqIndex];
		uint ToCol = 0;
		for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			char c = Seq[ColIndex];
			if (!Keeps[ColIndex])
				{
				asserta(isgap(c));
				continue;
				}
			Seq[ToCol++] = c;
			}
		asserta(ToCol == KeepCount);
		}
	M.m_uColCount = KeepCount;
	}

void TrimToRef(const MSA &Test, const MSA &Ref, MSA &Trimmed)
	{
	Trimmed.Clear();
	const uint TestSeqCount = Test.GetSeqCount();
	const uint TestColCount = Test.GetColCount();

	const uint RefSeqCount = Test.GetSeqCount();
	const uint RefColCount = Ref.GetColCount();

	map<string, uint> RefLabelToSeqIndex;
	map<string, uint> TestLabelToSeqIndex;
	vector<string> RefLabels;
	vector<string> TestLabels;
	Ref.GetLabelToSeqIndex(RefLabels, RefLabelToSeqIndex);
	Test.GetLabelToSeqIndex(TestLabels, TestLabelToSeqIndex);
	asserta(SIZE(RefLabels) == RefSeqCount);
	asserta(SIZE(TestLabels) == TestSeqCount);

	vector<string> Labels;
	vector<uint> RefSeqIndexes;
	for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
		{
		const string &Label = TestLabels[TestSeqIndex];
		map<string, uint>::const_iterator p =
		  RefLabelToSeqIndex.find(Label);
		if (p != RefLabelToSeqIndex.end())
			{
			uint RefSeqIndex = p->second;
			Labels.push_back(Label);
			RefSeqIndexes.push_back(RefSeqIndex);
			}
		}
	const uint TrimmedSeqCount = SIZE(RefSeqIndexes);

	Trimmed.SetSize(TrimmedSeqCount, TestColCount);
	for (uint TrimmedSeqIndex = 0; TrimmedSeqIndex < TrimmedSeqCount;
	  ++TrimmedSeqIndex)
		{
		uint RefSeqIndex = RefSeqIndexes[TrimmedSeqIndex];
		const string &Label = Labels[TrimmedSeqIndex];
		Trimmed.m_szNames[TrimmedSeqIndex] = mystrsave(Label.c_str());

		vector<bool> PosToUpper;
		const char *RefRow = Ref.m_szSeqs[RefSeqIndex];
		for (uint RefColIndex = 0; RefColIndex < RefColCount;
		  ++RefColIndex)
			{
			char c = RefRow[RefColIndex];
			if (!isgap(c))
				{
				bool Upper = isupper(c);
				PosToUpper.push_back(Upper);
				}
			}

		map<string, uint>::const_iterator p =
		  TestLabelToSeqIndex.find(Label);
		asserta(p != TestLabelToSeqIndex.end());
		uint TestSeqIndex = p->second;
		const char *TestRow = Test.m_szSeqs[TestSeqIndex];
		char *TrimmedRow = Trimmed.m_szSeqs[TrimmedSeqIndex];
		uint Pos = 0;
		for (uint TestColIndex = 0; TestColIndex < TestColCount;
		  ++TestColIndex)
			{
			char c = TestRow[TestColIndex];
			if (!isgap(c))
				{
				bool Upper = PosToUpper[Pos++];
				if (!Upper)
					c = '-';
				}
			TrimmedRow[TestColIndex] = c;
			}
		asserta(Pos == SIZE(PosToUpper));
		}
	DeleteAllGapColumns(Trimmed);
	}

void cmd_trimtoref()
	{
	MSA Test;
	MSA Ref;
	Test.FromFASTAFile(opt(trimtoref));
	Ref.FromFASTAFile_PreserveCase(opt(ref));

	MSA Trimmed;
	TrimToRef(Test, Ref, Trimmed);

	Trimmed.ToFASTAFile(opt(output));
	}
