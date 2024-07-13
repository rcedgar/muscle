#include "muscle.h"

// O(NL) computation of PREFAB Q score and Balibase TC score.
// Algorithm based on an idea due to Chuong (Tom) Do.
// Each position in the reference alignment is annotated with
// the column number C in the test alignment where the same 
// letter is found. A pair of identical Cs in the same reference
// column indicates a correctly aligned pair of letters.

void cmd_qscore()
	{
	const string TestFileName = opt(qscore);
	const string RefFileName = opt(ref);

	MSA msaTest;
	MSA msaRef;
	msaTest.FromFASTAFile(TestFileName);

	extern bool g_FASTA_Upper;
	bool SaveUpper = g_FASTA_Upper;
	g_FASTA_Upper = false;
	msaRef.FromFASTAFile(RefFileName);
	g_FASTA_Upper = SaveUpper;

	double Q = 0;
	double TC = 0;
	uint SeqDiffCount = 0;

	uint64 CorrectPairCount = 0;
	uint64 RefAlignedPairCount = 0;
	if (opt(verbose))
		{
		Log("RefCol  RefAln  NonGapped  TestAll  CorrCols  Ref\n");
		Log("------  ------  ---------  -------  --------  ---\n");
		//        6          9        7         8
		}

	const uint RefSeqCount = msaRef.GetSeqCount();
	const uint TestSeqCount = msaTest.GetSeqCount();

	const uint RefColCount = msaRef.GetColCount();
	const uint TestColCount = msaTest.GetColCount();

	map<string, uint> RefSeqNameToIndex;
	map<string, uint> RefSeqToIndex;
	vector<uint> RefToTestSeqIndex(RefSeqCount);

	for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const string SeqName = msaRef.GetSeqName(RefSeqIndex);
		string USeq;
		msaRef.GetUngappedSeqStr(RefSeqIndex, USeq);
		RefToTestSeqIndex[RefSeqIndex] = UINT_MAX;
		RefSeqNameToIndex[SeqName] = RefSeqIndex;
		RefSeqToIndex[USeq] = RefSeqIndex;
		}

	uint FoundCount = 0;
	for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
		{
		const string SeqName = msaTest.GetSeqName(TestSeqIndex);
		if (opt(bysequence))
			{
			string USeq;
			msaTest.GetUngappedSeqStr(TestSeqIndex, USeq);
			map<string, uint>::const_iterator p =
			  RefSeqToIndex.find(USeq);
			if (p != RefSeqToIndex.end())
				{
				uint RefSeqIndex = p->second;
				if (RefSeqIndex == UINT_MAX)
					Die("UINT_MAX");
				RefToTestSeqIndex[RefSeqIndex] = TestSeqIndex;
				++FoundCount;
				}
			}
		else
			{
			map<string, uint>::const_iterator p =
			  RefSeqNameToIndex.find(SeqName);
			if (p != RefSeqNameToIndex.end())
				{
				uint RefSeqIndex = p->second;
				if (RefSeqIndex == UINT_MAX)
					Die("UINT_MAX");
				RefToTestSeqIndex[RefSeqIndex] = TestSeqIndex;
				++FoundCount;
				}
			}
		}
	if (opt(bysequence))
		{
		if (FoundCount < 2)
			Die("%d ref seqs found in test MSA", FoundCount);
		if (FoundCount < RefSeqCount)
			Warning("%u / %u ref seqs found in test MSA",
				FoundCount, RefSeqCount);
		}
	else
		{
		if (FoundCount < 2)
			Die("%d ref labels found in test MSA", FoundCount);
		if (FoundCount < RefSeqCount)
			Warning("%u / %u ref labels found in test MSA",
				FoundCount, RefSeqCount);
		}

// TestColIndex[i] is the one-based (not zero-based!) test column index
// of the letter found in the current column of the reference alignment
// (or the most recent letter if the reference column is gapped, or zero
// if no letter has yet been found). Here, seq index i is for msaRef.
	vector<uint> TestColIndex(TestSeqCount, 0);

// TestColIndexCount[i] is the number of times that a letter from test
// column i (one-based!) appears in the current reference column.
	vector<uint> TestColIndexCount(TestColCount+1, 0);

// TestColIndexes[i] is the column index in the test alignment of
// the i'th non-gapped position in the current reference column.
	vector<uint> TestColIndexes;

	uint RefAlignedColCount = 0;
	uint CorrectColCount = 0;
	if (opt(verbose))
		{
		Log("RefCol  RefAln  NonGapped  TestAll  CorrCols  Ref\n");
		Log("------  ------  ---------  -------  --------  ---\n");
		//        6          9        7         8
		}

	for (uint RefColIndex = 0; RefColIndex < RefColCount; RefColIndex++)
		{
		TestColIndexes.clear();
		TestColIndexes.reserve(RefSeqCount);

	// NonGappedCount is the number of non-gapped positions in the current
	// reference column.
		uint NonGappedCount = 0;
		uint FirstTestColIndex = UINT_MAX;
		bool RefColIsAligned = false;
		bool TestColAllCorrect = true;
		bool TestAllAligned = true;
		for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
			{
			uint TestSeqIndex = RefToTestSeqIndex[RefSeqIndex];
			if (TestSeqIndex == UINT_MAX)
				continue;

			char cRef = msaRef.GetChar(RefSeqIndex, RefColIndex);
			if (!isgap(cRef))
				{
				char cTest = 0;
				uint Col = TestColIndex[TestSeqIndex];
				do
					cTest = msaTest.GetChar(TestSeqIndex, Col++);
				while (isgap(cTest));
				if (toupper(cRef) != toupper(cTest))
					{
					++SeqDiffCount;
					Warning("Test seq %u (%s) differs from ref seq %u (%s), ref col %u=%c, test=%c",
						TestSeqIndex,
						msaTest.GetSeqName(TestSeqIndex),
						RefSeqIndex,
						msaRef.GetSeqName(RefSeqIndex),
						RefColIndex,
						cRef,
						cTest);
					}
				if (isalpha(cRef) && isupper(cRef))
					{
					RefColIsAligned = true;
					++NonGappedCount;
					if (isupper(cTest))
						{
						TestColIndexes.push_back(Col);
						++(TestColIndexCount[Col]);
						if (FirstTestColIndex == UINT_MAX)
							FirstTestColIndex = Col;
						else
							{
							if (FirstTestColIndex != Col)
								TestColAllCorrect = false;
							}
						}
					else
						TestAllAligned = false;
					}
				else
					{
					if (RefColIsAligned)
						{
						Log("\n");
						Log("Ref col: ");
						for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
							Log("%c", msaRef.GetChar(RefSeqIndex, RefColIndex));
						Log("\n");
						Die("Ref col %u has both upper- and lower-case letters",
						  RefColIndex);
						}
					}
				TestColIndex[TestSeqIndex] = Col;
				}
			}

		if (RefColIsAligned && NonGappedCount > 1)
			{
			++RefAlignedColCount;
			if (TestColAllCorrect && TestAllAligned)
				++CorrectColCount;
			}

		uint ColPairCount = 0;
		for (vector<uint>::const_iterator p = TestColIndexes.begin(); p != TestColIndexes.end(); ++p)
			{
			uint Col = *p;
			uint Count = TestColIndexCount[Col];
			if (Count > 0)
				ColPairCount += Count*(Count - 1)/2;
			TestColIndexCount[Col] = 0;
			}

		CorrectPairCount += ColPairCount;
		RefAlignedPairCount += NonGappedCount*(NonGappedCount - 1)/2;

		if (opt(verbose))
			{
			Log("%6u  %6c  %9u  %7c  %8u  ",
			  RefColIndex, 
			  RefColIsAligned ? 'T' : 'F',
			  NonGappedCount,
			  TestColAllCorrect ? 'T' : 'F',
			  CorrectColCount);
			for (uint RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
				{
				uint TestSeqIndex = RefToTestSeqIndex[RefSeqIndex];
				if (TestSeqIndex == UINT_MAX)
					continue;

				char cRef = msaRef.GetChar(RefSeqIndex, RefColIndex);
				Log("%c", cRef);
				}
			Log("\n");
			}
		}

	if (RefAlignedPairCount == 0)
		Q = 0;
	else
		Q = (double) CorrectPairCount / (double) RefAlignedPairCount;

	if (RefAlignedColCount == 0)
		{
		Warning("reference alignment %s has no aligned (upper-case) columns\n",
		  RefFileName.c_str());
		TC = 0;
		}
	else
		TC = (double) CorrectColCount / (double) RefAlignedColCount;

	if (opt(verbose))
		{
		Log("        ------                      --------\n");
		Log("%6.6s  %6u  %9.9s  %7.7s  %8u\n",
		  "", RefAlignedColCount, "", "", CorrectColCount);
		Log("\n");
		Log("CorrectPairCount     %u\n", CorrectPairCount);
		Log("RefAlignedPairCount  %u\n", RefAlignedPairCount);
		Log("CorrectColCount      %u\n", CorrectColCount);
		Log("RefAlignedColCount   %u\n", RefAlignedColCount);
		Log("Q                    %.4f\n", Q);
		Log("TC                   %.4f\n", TC);
		}

	if (SeqDiffCount > 0)
		Warning("%u seq diffs ignored", SeqDiffCount);

	ProgressLog("%s Q=%.3g, TC=%.3g\n", TestFileName.c_str(), Q, TC);
	}
