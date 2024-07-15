#include "muscle.h"
#include "mega.h"

float *CalcPost(uint GSIX, uint GSIY)
	{
	uint LX = GetSeqLengthByGSI(GSIX);
	uint LY = GetSeqLengthByGSI(GSIY);
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);

	float *Fwd = AllocFB(LX, LY);
	float *Bwd = AllocFB(LX, LY);

	if (Mega::m_Loaded)
		{
		const vector<vector<byte> > &ProfileX = *Mega::GetProfileByGSI(GSIX);
		const vector<vector<byte> > &ProfileY = *Mega::GetProfileByGSI(GSIY);
		asserta(SIZE(ProfileX) == LX);
		asserta(SIZE(ProfileY) == LY);
		Mega::CalcFwdFlat_mega(ProfileX, ProfileY, Fwd);
		Mega::CalcBwdFlat_mega(ProfileX, ProfileY, Bwd);
		}
	else
		{
		const byte *X = GetByteSeqByGSI(GSIX);
		const byte *Y = GetByteSeqByGSI(GSIY);
		CalcFwdFlat(X, LX, Y, LY, Fwd);
		CalcBwdFlat(X, LX, Y, LY, Bwd);
		}

	float *Post = AllocPost(LX, LY);
	CalcPostFlat(Fwd, Bwd, LX, LY, Post);
	myfree(Fwd);
	myfree(Bwd);
	return Post;
	}
