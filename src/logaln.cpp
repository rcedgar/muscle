#include "muscle.h"
#include "pathinfo.h"

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
		if (c == 'B' || c == 'M')
			{
			RowX += XSeq[XPos];
			RowY += YSeq[YPos];
			++YPos;
			++XPos;
			}
		else if (c == 'X' || c == 'D')
			{
			RowX += XSeq[XPos];
			RowY += '-';
			++XPos;
			}
		else if (c == 'Y' || c == 'I')
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

void MakeAlnRows(const byte *XSeq, uint LX,
  const byte *YSeq, uint LY, const PathInfo &PI,
  string &RowX, string &RowY)
	{
	RowX.clear();
	RowY.clear();

	string PathXY;
	PI.GetPathStr(PathXY);
	const uint ColCount = SIZE(PathXY);
	uint XPos = 0;
	uint YPos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = PathXY[Col];
		if (c == 'B' || c == 'M')
			{
			RowX += XSeq[XPos];
			RowY += YSeq[YPos];
			++YPos;
			++XPos;
			}
		else if (c == 'X' || c == 'D')
			{
			RowX += XSeq[XPos];
			RowY += '-';
			++XPos;
			}
		else if (c == 'Y' || c == 'I')
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

void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI)
	{
	string sX;
	string sY;
	for (uint i = 0; i < LX; ++i)
		sX += X[i];
	for (uint i = 0; i < LY; ++i)
		sY += Y[i];
	string PathXY;
	PI.GetPathStr(PathXY);
	LogAln(sX, sY, PathXY);
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

void WriteAnnotRow(FILE *f, const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi)
	{
	fprintf(f, "%5.5s ", "");
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M')
			{
			byte a = A[i++];
			byte b = B[j++];
			if (toupper(a) == toupper(b))
				fprintf(f, "|");
			else
				fprintf(f, " ");
			}
		else
			{
			if (c == 'D')
				++i;
			else if (c == 'I')
				++j;
			else
				asserta(false);
			fprintf(f, " ");
			}
		}
	fprintf(f, "\n");
	}

void WriteBRow(FILE *f, const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi, const string &LabelB)
	{
	fprintf(f, "%5u ", j+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'I')
			fprintf(f, "%c", B[j++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u  %s\n", j, LabelB.c_str());
	}

void WriteARow(FILE *f, const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi, const string &LabelA)
	{
	fprintf(f, "%5u ", i+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'D')
			fprintf(f, "%c", A[i++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u  %s\n", i, LabelA.c_str());
	}

void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path)
	{
	const unsigned BLOCK_SIZE = 80;
	unsigned ALo, BLo, ColLo, ColHi;
	ALo = 0;
	BLo = 0;
	ColLo = 0;
	ColHi = (unsigned) strlen(Path) - 1;

	asserta(ColHi >= ColLo);

	unsigned i = ALo;
	unsigned j = BLo;
	unsigned ColFrom = ColLo;
	for (;;)
		{
		if (ColFrom > ColHi)
			break;
		unsigned ColTo = ColFrom + BLOCK_SIZE - 1;
		if (ColTo > ColHi)
			ColTo = ColHi;

		unsigned i0 = i;
		unsigned j0 = j;
		WriteARow(f, A, Path, i, ColFrom, ColTo, "");
		WriteAnnotRow(f, A, B, Path, i0, j0, ColFrom, ColTo);
		WriteBRow(f, B, Path, j, ColFrom, ColTo, "");
		fprintf(f, "\n");

		ColFrom += BLOCK_SIZE;
		}
	}
