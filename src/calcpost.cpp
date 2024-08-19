#include "muscle.h"
#include "mega.h"

float *CalcPost(const string &LabelX, const string &LabelY)
	{
	uint LX = GetSeqLengthByGlobalLabel(LabelX);
	uint LY = GetSeqLengthByGlobalLabel(LabelY);
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);

	float *Fwd = AllocFB(LX, LY);
	float *Bwd = AllocFB(LX, LY);

	if (Mega::m_Loaded)
		{
		const vector<vector<byte> > &ProfileX = *Mega::GetProfileByLabel(LabelX);
		const vector<vector<byte> > &ProfileY = *Mega::GetProfileByLabel(LabelY);
		asserta(SIZE(ProfileX) == LX);
		asserta(SIZE(ProfileY) == LY);
		Mega::CalcFwdFlat_mega(ProfileX, ProfileY, Fwd);
		Mega::CalcBwdFlat_mega(ProfileX, ProfileY, Bwd);
		}
	else
		{
		const byte *X = GetGlobalByteSeqByLabel(LabelX);
		const byte *Y = GetGlobalByteSeqByLabel(LabelY);
		CalcFwdFlat(X, LX, Y, LY, Fwd);
		CalcBwdFlat(X, LX, Y, LY, Bwd);
		}

	float *Post = AllocPost(LX, LY);
	CalcPostFlat(Fwd, Bwd, LX, LY, Post);
	myfree(Fwd);
	myfree(Bwd);
	return Post;
	}

//float *CalcPost(const string &Label1, const string &Label2)
//	{
//	uint GSI1 = GetGSIByLabel(Label1);
//	uint GSI2 = GetGSIByLabel(Label2);
//	float *Post = CalcPost(GSI1, GSI2);
//	return Post;
//	}
