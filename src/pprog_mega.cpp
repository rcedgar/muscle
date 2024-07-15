#include "muscle.h"
#include "pprog_mega.h"

//void CalcFwdFlat_mega(const Mega &M,
//  const vector<vector<byte> > &ProfileX,
//  const vector<vector<byte> > &ProfileY, float *Flat)

void PProg_mega::CalcFwdFlat_PProg(uint GSI1, uint L1, 
  uint GSI2, uint L2, float *Flat)
	{
	const vector<vector<byte> > &Profile1 = *Mega::GetProfileByGSI(GSI1);
	const vector<vector<byte> > &Profile2 = *Mega::GetProfileByGSI(GSI2);
	asserta(SIZE(Profile1) == L1);
	asserta(SIZE(Profile2) == L2);
	Mega::CalcFwdFlat_mega(Profile1, Profile2, Flat);
	}

void PProg_mega::CalcBwdFlat_PProg(uint GSI1, uint L1, 
  uint GSI2, uint L2, float *Flat)
	{
	const vector<vector<byte> > &Profile1 = *Mega::GetProfileByGSI(GSI1);
	const vector<vector<byte> > &Profile2 = *Mega::GetProfileByGSI(GSI2);
	asserta(SIZE(Profile1) == L1);
	asserta(SIZE(Profile2) == L2);
	Mega::CalcBwdFlat_mega(Profile1, Profile2, Flat);
	}
