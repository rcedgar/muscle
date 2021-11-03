#include "muscle.h"

void MakeAlnRows(const string &X, const string &Y,
  const string &PathXY, string &RowX, string &RowY)
	{
	RowX.clear();
	RowY.clear();

	const byte *XSeq = (const byte *) X.c_str();
	const byte *YSeq = (const byte *) Y.c_str();
	const uint ColCount = SIZE(PathXY);
	const uint LX = SIZE(X);
	const uint LY = SIZE(Y);
	uint XPos = 0;
	uint YPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = PathXY[Col];
		if (c == 'B')
			{
			RowX += XSeq[XPos];
			RowY += YSeq[YPos];
			++YPos;
			++XPos;
			}
		else if (c == 'X')
			{
			RowX += XSeq[XPos];
			RowY += '-';
			++XPos;
			}
		else if (c == 'Y')
			{
			RowY += YSeq[YPos];
			RowX += '-';
			++YPos;
			}
		else
			asserta(false);
		}
	asserta(XPos == LX && YPos == LY);
	}

void MakeAlnRows(const Sequence &X, const Sequence &Y,
  const string &PathXY, string &RowX, string &RowY)
	{
	RowX.clear();
	RowY.clear();

	const byte *XSeq = X.GetBytePtr();
	const byte *YSeq = Y.GetBytePtr();
	const uint ColCount = SIZE(PathXY);
	const uint LX = X.GetLength();
	const uint LY = Y.GetLength();
	uint XPos = 0;
	uint YPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = PathXY[Col];
		if (c == 'B')
			{
			RowX += XSeq[XPos];
			RowY += YSeq[YPos];
			++YPos;
			++XPos;
			}
		else if (c == 'X')
			{
			RowX += XSeq[XPos];
			RowY += '-';
			++XPos;
			}
		else if (c == 'Y')
			{
			RowY += YSeq[YPos];
			RowX += '-';
			++YPos;
			}
		else
			asserta(false);
		}
	asserta(XPos == LX && YPos == LY);
	}

void LogAln(const string &X, const string &Y, const string &PathXY)
	{
	string RowX;
	string RowY;
	MakeAlnRows(X, Y, PathXY, RowX, RowY);
	Log("\n");
	Log("%s\n", RowX.c_str());
	Log("%s\n", RowY.c_str());
	}

void LogAln(const Sequence &X, const Sequence &Y, const string &PathXY)
	{
	string RowX;
	string RowY;
	MakeAlnRows(X, Y, PathXY, RowX, RowY);
	const string &LabelX = X.GetLabel();
	const string &LabelY = Y.GetLabel();
	Log("\n");
	Log("%10.10s  %s\n", LabelX.c_str(), RowX.c_str());
	Log("%10.10s  %s\n", LabelY.c_str(), RowY.c_str());
	}

void PathToColVecs(const string &PathXY,
  vector<uint> &PosToColX, vector<uint> &PosToColY,
  vector<uint> &ColToPosX, vector<uint> &ColToPosY)
	{
	PosToColX.clear();
	PosToColY.clear();
	ColToPosX.clear();
	PosToColY.clear();

	const uint ColCount = SIZE(PathXY);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = PathXY[Col];
		if (c == 'B')
			{
			uint PosX = SIZE(PosToColX);
			uint PosY = SIZE(PosToColY);

			ColToPosX.push_back(PosX);
			ColToPosY.push_back(PosY);

			PosToColX.push_back(Col);
			PosToColY.push_back(Col);
			}
		else if (c == 'X')
			{
			uint PosX = SIZE(PosToColX);
			ColToPosX.push_back(PosX);
			ColToPosY.push_back(UINT_MAX);

			PosToColX.push_back(Col);
			}
		else if (c == 'Y')
			{
			uint PosY = SIZE(PosToColY);
			ColToPosY.push_back(PosY);
			ColToPosX.push_back(UINT_MAX);

			PosToColY.push_back(Col);
			}
		else
			asserta(false);
		}

// Validate
	{
	asserta(SIZE(ColToPosX) == ColCount);
	asserta(SIZE(ColToPosY) == ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint PosX = ColToPosX[Col];
		uint PosY = ColToPosY[Col];
		if (PosX != UINT_MAX)
			{
			asserta(PosX < SIZE(PosToColX));
			asserta(PosToColX[PosX] == Col);
			}
		if (PosY != UINT_MAX)
			{
			asserta(PosY < SIZE(PosToColY));
			asserta(PosToColY[PosY] == Col);
			}
		}
	}
	}
