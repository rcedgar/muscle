#include "muscle.h"
#include "mpcflat.h"

void RelaxFlat_XZ_ZY(const MySparseMx &XZ, const MySparseMx &ZY,
  float WeightZ, float *Post)
	{
	const uint LX = XZ.GetLX();
	const uint LZ = XZ.GetLY();
	const uint LY = ZY.GetLY();
	asserta(ZY.GetLX() == LZ);

	for (uint PosX = 0; PosX < LX; ++PosX)
		{
		uint Offset_XZ = XZ.GetOffset(PosX);
		uint Size_XZ = XZ.GetSize(PosX);
		for (uint k = 0; k < Size_XZ; ++k)
			{
			float P_XZ = XZ.GetProb_Offset(Offset_XZ + k);
			uint PosZ = XZ.GetCol_Offset(Offset_XZ + k);

			uint Offset_ZY = ZY.GetOffset(PosZ);
			uint Size_ZY = ZY.GetSize(PosZ);
			for (uint m = 0; m < Size_ZY; ++m)
				{
				float P_ZY = ZY.GetProb_Offset(Offset_ZY + m);
				uint PosY = ZY.GetCol_Offset(Offset_ZY + m);
				Post[PosX*LY + PosY] += WeightZ*P_XZ*P_ZY;
				}
			}
		}
	}

void RelaxFlat_ZX_ZY(const MySparseMx &ZX, const MySparseMx &ZY,
  float WeightZ, float *Post)
	{
	const uint LZ = ZX.GetLX();
	const uint LX = ZX.GetLY();
	const uint LY = ZY.GetLY();
	asserta(ZY.GetLX() == LZ);

	for (uint PosZ = 0; PosZ < LZ; ++PosZ)
		{
		uint Offset_ZX = ZX.GetOffset(PosZ);
		uint Size_ZX = ZX.GetSize(PosZ);
		for (uint k = 0; k < Size_ZX; ++k)
			{
			float P_ZX = ZX.GetProb_Offset(Offset_ZX + k);
			uint PosX = ZX.GetCol_Offset(Offset_ZX + k);

			uint Offset_ZY = ZY.GetOffset(PosZ);
			uint Size_ZY = ZY.GetSize(PosZ);
			for (uint m = 0; m < Size_ZY; ++m)
				{
				float P_ZY = ZY.GetProb_Offset(Offset_ZY + m);
				uint PosY = ZY.GetCol_Offset(Offset_ZY + m);
				Post[PosX*LY + PosY] += WeightZ*P_ZX*P_ZY;
				}
			}
		}
	}

void RelaxFlat_XZ_YZ(const MySparseMx &XZ, const MySparseMx &YZ,
  float WeightZ, float *Post)
	{
	const uint LX = XZ.GetLX();
	const uint LZ = XZ.GetLY();
	const uint LY = YZ.GetLX();
	asserta(YZ.GetLY() == LZ);

	vector<uint> PosZToLoPosY;
	vector<uint> PosZToHiPosY;
	YZ.GetColToRowLoHi(PosZToLoPosY, PosZToHiPosY);

	for (uint PosX = 0; PosX < LX; ++PosX)
		{
		uint Offset_XZ = XZ.GetOffset(PosX);
		uint Size_XZ = XZ.GetSize(PosX);
		for (uint k = 0; k < Size_XZ; ++k)
			{
			float P_XZ = XZ.GetProb_Offset(Offset_XZ + k);
			uint PosZ = XZ.GetCol_Offset(Offset_XZ + k);

			uint LoPosY = PosZToLoPosY[PosZ];
			uint HiPosY = PosZToHiPosY[PosZ];
			if (LoPosY == UINT_MAX)
				continue;
			for (uint PosY = LoPosY; PosY <= HiPosY; ++PosY)
				{
				float P_YZ = YZ.GetProb(PosY, PosZ);
				Post[PosX*LY + PosY] += WeightZ*P_XZ*P_YZ;
				}
			}
		}
	}
