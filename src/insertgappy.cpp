#include "myutils.h"

void InsertGappyPositions(const vector<char> &OccPath,
  uint FullColCount1, uint FullColCount2,
  vector<uint> &OccCols1, vector<uint> &OccCols2,
  vector<char> &FullPath)
	{
	FullPath.clear();

	const uint CombinedOccColCount = SIZE(OccPath);
	const uint OccColCount1 = SIZE(OccCols1);
	const uint OccColCount2 = SIZE(OccCols2);

	uint CombinedOccCol = 0;
	uint OccPos1 = 0;
	uint OccPos2 = 0;
	uint FullPos1 = 0;
	uint FullPos2 = 0;
	for (uint CombinedOccCol = 0; CombinedOccCol < CombinedOccColCount; ++CombinedOccCol)
		{
		char OccPathOp = OccPath[CombinedOccCol];
		bool Add1 = false;
		bool Add2 = false;
		switch (OccPathOp)
			{
		case 'X':
			Add1 = true;
			Add2 = false;
			break;
		case 'Y':
			Add1 = false;
			Add2 = true;
			break;
		case 'B':
			Add1 = true;
			Add2 = true;
			break;
			}
		uint InsertCount1 = 0;
		uint InsertCount2 = 0;
		if (Add1 && OccPos1 + 1 < OccColCount1)
			{
			uint NewFullPos1 = OccCols1[OccPos1];
			asserta(NewFullPos1 >= FullPos1);
			InsertCount1 = NewFullPos1 - FullPos1;
			}
		if (Add2 && OccPos2 + 1 < OccColCount2)
			{
			uint NewFullPos2 = OccCols2[OccPos2];
			asserta(NewFullPos2 >= FullPos2);
			InsertCount2 = NewFullPos2 - FullPos2;
			}
		uint MinInsertCount = min(InsertCount1, InsertCount2);
		for (uint i = 0; i < MinInsertCount; ++i)
			{
			FullPath.push_back('B');
			++FullPos1;
			++FullPos2;
			}
		for (uint i = MinInsertCount; i < InsertCount1; ++i)
			{
			if (!FullPath.empty() && FullPath.back() == 'Y')
				FullPath.back() = 'B';
			else
				FullPath.push_back('X');
			++FullPos1;
			}
		for (uint i = MinInsertCount; i < InsertCount2; ++i)
			{
			if (!FullPath.empty() && FullPath.back() == 'X')
				FullPath.back() = 'B';
			else
				FullPath.push_back('Y');
			++FullPos2;
			}

		FullPath.push_back(OccPathOp);
		if (Add1)
			{
			++FullPos1;
			++OccPos1;
			}
		if (Add2)
			{
			++FullPos2;
			++OccPos2;
			}
		}

	uint InsertCount1 = 0;
	uint InsertCount2 = 0;
	if (FullPos1 < FullColCount1)
		InsertCount1 = FullColCount1 - FullPos1;
	if (FullPos2 < FullColCount2)
		InsertCount2 = FullColCount2 - FullPos2;
	uint MinInsertCount = min(InsertCount1, InsertCount2);
	for (uint i = 0; i < MinInsertCount; ++i)
		{
		FullPath.push_back('B');
		++FullPos1;
		++FullPos2;
		}
	for (uint i = MinInsertCount; i < InsertCount1; ++i)
		{
		if (!FullPath.empty() && FullPath.back() == 'Y')
			FullPath.back() = 'B';
		else
			FullPath.push_back('X');
		++FullPos1;
		}
	for (uint i = MinInsertCount; i < InsertCount2; ++i)
		{
		if (!FullPath.empty() && FullPath.back() == 'X')
			FullPath.back() = 'B';
		else
			FullPath.push_back('Y');
		++FullPos2;
		}

	asserta(OccPos1 == OccColCount1);
	asserta(OccPos2 == OccColCount2);

	asserta(FullPos1 == FullColCount1);
	asserta(FullPos2 == FullColCount2);
	}
