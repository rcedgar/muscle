#include "myutils.h"
#include "msa.h"

// MergeMap[i] is column in X2 and Y2 where the i'th
//  entry in ColsX and ColsY is placed in Path.
void AlignMSAsByCols(const MSA &X, const MSA &Y,
  const vector<uint> &ColsX, const vector<uint> &ColsY,
  string &Path, vector<uint> &MergeMap, MSA &X2, MSA &Y2)
	{
	X2.Clear();
	Y2.Clear();
	MergeMap.clear();
	const uint NC = SIZE(ColsX);
	asserta(SIZE(ColsY) == NC);
	asserta(NC > 0);
	for (uint i = 1; i < NC; ++i)
		{
		asserta(ColsX[i] > ColsX[i-1]);
		asserta(ColsY[i] > ColsY[i-1]);
		}

	const uint ColCountX = X.GetColCount();
	const uint ColCountY = Y.GetColCount();
	uint ColX = 0;
	uint ColY = 0;
	uint FirstColX = ColsX[0];
	uint FirstColY = ColsY[0];
	if (FirstColX < FirstColY)
		{
		for (uint i = 0; i < FirstColY - FirstColX; ++i)
			{
			Path += 'y';
			++ColY;
			}
		}
	else if (FirstColY < FirstColX)
		{
		for (uint i = 0; i < FirstColX - FirstColY; ++i)
			{
			Path += 'x';
			++ColX;
			}
		}
	while (ColX < FirstColX)
		{
		Path += 'm';
		++ColX;
		++ColY;
		}
	asserta(ColX == FirstColX && ColY == FirstColY);
	MergeMap.push_back(SIZE(Path));
	Path += 'M';
	for (uint i = 1; i < NC; ++i)
		{
		uint NextColX = ColsX[i];
		uint NextColY = ColsY[i];
		asserta(ColX < NextColX);
		asserta(ColY < NextColY);
		uint nx = NextColX - ColX;
		uint ny = NextColY - ColY;
		if (nx < ny)
			{
			for (uint j = 0; j < ny-nx; ++j)
				{
				Path += 'y';
				++ColY;
				}
			}
		else if (ny < nx)
			{
			for (uint j = 0; j < nx-ny; ++j)
				{
				Path += 'x';
				++ColX;
				}
			}
		while (ColX + 1 < NextColX)
			{
			Path += 'm';
			++ColX;
			++ColY;
			}
		asserta(ColX + 1 == NextColX);
		asserta(ColY + 1 == NextColY);
		MergeMap.push_back(SIZE(Path));
		Path += 'M';
		ColX = NextColX;
		ColY = NextColY;
		}
	++ColX;
	++ColY;
	while (ColX < ColCountX && ColY < ColCountY)
		{
		Path += 'm';
		++ColX;
		++ColY;
		}
	while (ColX < ColCountX)
		{
		Path += 'x';
		++ColX;
		}
	while (ColY < ColCountY)
		{
		Path += 'y';
		++ColY;
		}

	const uint PL = SIZE(Path);
	uint nx = 0;
	uint ny = 0;
	for (uint i = 0; i < PL; ++i)
		{
		char c = Path[i];
		switch (c)
			{
		case 'm': case 'M': ++nx; ++ny; break;
		case 'x': ++nx; break;
		case 'y': ++ny; break;
		default: asserta(false);
			}
		}
	asserta(nx == ColCountX);
	asserta(ny == ColCountY);

	const uint SeqCountX = X.GetSeqCount();
	//uint Sz = LastM - FirstM + 1;
	X2.SetSize(SeqCountX, PL);
	for (uint i = 0; i < SeqCountX; ++i)
		{
		X2.m_szNames[i] = mystrsave(X.m_szNames[i]);
		uint Col = 0;
		for (uint j = 0; j < PL; ++j)
			{
			char c = Path[j];
			switch (c)
				{
			case 'm':
			case 'M':
			case 'x':
				{
				char c = X.GetChar(i, Col);
				X2.SetChar(i, j, c);
				++Col;
				break;
				}

			case 'y':
				{
				X2.SetChar(i, j, '.');
				break;
				}

			default:	asserta(false);
				}
			}
		}

	const uint SeqCountY = Y.GetSeqCount();
	Y2.SetSize(SeqCountY, PL);
	for (uint i = 0; i < SeqCountY; ++i)
		{
		Y2.m_szNames[i] = mystrsave(Y.m_szNames[i]);
		uint Col = 0;
		for (uint j = 0; j < PL; ++j)
			{
			char c = Path[j];
			switch (c)
				{
			case 'm':
			case 'M':
			case 'y':
				{
				char c = Y.GetChar(i, Col);
				Y2.SetChar(i, j, c);
				++Col;
				break;
				}

			case 'x':
				{
				Y2.SetChar(i, j, '.');
				break;
				}

			default:	asserta(false);
				}
			}
		}
	//X2.LogMe();
	//Y2.LogMe();
	}
