#pragma once

struct DiagBox;

void GetDiagBox(uint LA, uint LB, uint DiagLo, uint DiagHi, DiagBox &Box);
void GetDiagRange(uint LA, uint LB, uint d,
  uint &mini, uint &minj, uint &maxi, uint &maxj);
void GetDiagLoHi(uint LA, uint LB, const char *Path,
  uint &dlo, uint &dhi);

struct DiagBox
	{
	DiagBox()
		{
		}

	DiagBox(uint LA_, uint LB_, uint DiagLo, uint DiagHi)
		{
		//GetDiagBox(LA, LB, DiagLo, DiagHi, *this);
		//Validate();
		Init(LA_, LB_, DiagLo, DiagHi);
		}

	void Init(uint LA_, uint LB_, uint DiagLo, uint DiagHi)
		{
		GetDiagBox(LA_, LB_, DiagLo, DiagHi, *this);
		Validate();
		}

	uint LA;
	uint LB;

	uint dlo;
	uint dhi;

	uint dlo_mini;
	uint dlo_minj;

	uint dlo_maxi;
	uint dlo_maxj;

	uint dhi_mini;
	uint dhi_minj;

	uint dhi_maxi;
	uint dhi_maxj;

	uint GetDiag(uint i, uint j) const
		{
		return LA - i + j;
		}

// i, j are positions 0..LA-1, 0..LB-1.
	bool InBox(uint i, uint j) const
		{
		uint d = GetDiag(i, j);
		return d >= dlo && d <= dhi;
		}

/***
i, j are 0-based prefix lengths 0..LA, 0..LB.

A full path is in the box iff all match pairs are in the box.

A partial path that aligns a prefix of A to a prefix of B as
in D.P.) is in the box iff it is is the prefix of at least
one full path that is in the box.

A D.P. matrix entry X[i][j] is in the box iff there is at
least one full path aligning the first i letters of A and
the first j letters of B ending in a column of type X, i.e.
if there exists a partial path in the box that ends in X.

Assume terminals appear in all paths, and DI/ID forbidden.

Intuitively seems that by these definitions D is in box iff
DM or MD is in box, I is in box iff IM or MI is in box.
Don't have proof..
***/
	bool InBoxDPM(uint i, uint j) const
		{
	// Special case for M[0][0]
		if (i == 0 && j == 0)
			return true;
		if (i == 0 || j == 0)
			return false;
		uint d = GetDiag(i-1, j-1);
		return d >= dlo && d <= dhi;
		}

	bool InBoxDPD(uint i, uint j) const
		{
		bool MD = i == 0 ? false : InBoxDPM(i-1, j);
		bool DM = (i == LA || j == LB) ? false : InBoxDPM(i+1, j+1);
		return MD || DM;
		}

	bool InBoxDPI(uint i, uint j) const
		{
		bool MI = j == 0 ? false : InBoxDPM(i, j-1);
		bool IM = (i == LA || j == LB) ? false : InBoxDPM(i+1, j+1);
		return MI || IM;
		}

	// d = LA - i + j = 1 .. LA+LB-1
	void Validate() const
		{
		asserta(dlo <= dhi);
		asserta(dlo >= GetDiag(LA-1, 0));
		asserta(dhi <= GetDiag(0, LB-1));

		asserta(GetDiag(dlo_mini, dlo_minj) == dlo);
		asserta(GetDiag(dlo_maxi, dlo_maxj) == dlo);
		asserta(GetDiag(dhi_mini, dhi_minj) == dhi);
		asserta(GetDiag(dhi_maxi, dhi_maxj) == dhi);

		asserta(dlo_mini >= dhi_mini);
		asserta(dlo_minj <= dhi_minj);
		asserta(dlo_maxi >= dhi_maxi);
		asserta(dlo_maxj <= dhi_maxj);
		}

	uint GetMini() const
		{
		return dhi_mini;
		}

	uint GetMaxi() const
		{
		return dlo_maxi;
		}

	uint GetMinj() const
		{
		return dlo_minj;
		}

	uint GetMaxj() const
		{
		return dhi_maxj;
		}
/***
	i = 0..LA-1
	j = 0..LB-1
	d = LA - i + j = 1 .. LA+LB-1
	j = d - LA + i
	i = LA - d + j
***/
	void GetRange_j(uint i, uint &Startj, uint &Endj) const
		{
	// j = d - LA + i
		if (dlo + i >= LA)
			Startj = dlo + i - LA;
		else
			Startj = 0;

		if (Startj >= LB)
			Startj = LB - 1;

		if (dhi + i + 1 >= LA)
			Endj = dhi + i + 1 - LA;
		else
			Endj = 0;

		if (Endj >= LB)
			Endj = LB - 1;

		asserta(Endj >= Startj);
		asserta(Startj < LB);
		}

	void LogMe() const
		{
		Log("LA=%u LB=%d dlo(%u): (%u,%u)-(%u,%u) dhi(%u): (%u,%u)-(%u,%u) i=[%u-%u] j=[%u-%u]\n",
		  LA, LB,
		  dlo,
		  dlo_mini, dlo_minj,
		  dlo_maxi, dlo_maxj,
		  dhi,
		  dhi_mini, dhi_minj,
		  dhi_maxi, dhi_maxj,
		  GetMini(), GetMaxi(),
		  GetMinj(), GetMaxj());
		}
	};
