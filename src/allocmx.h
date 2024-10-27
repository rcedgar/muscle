#pragma once

template <typename T> void AllocMx(vector<vector<T> > &Mx,
  uint LA, uint LB, T InitialValue)
	{
	Mx.clear();
	Mx.resize(LA);
	for (uint i = 0; i < LA; ++i)
		Mx[i].resize(LB, InitialValue);
	}

