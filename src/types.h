#pragma once

class MSA;
class Seq;
class TextFile;
class Tree;
class SeqVect;

enum LINKAGE
	{
	LINKAGE_Undefined,
	LINKAGE_Min,
	LINKAGE_Max,
	LINKAGE_Avg,
	LINKAGE_Biased,
	};

typedef float t_ByteVec[256];
typedef float t_ByteMx[256][256];

typedef float SCOREMATRIX[32][32];
typedef SCOREMATRIX *PTR_SCOREMATRIX;

typedef float Mx2020[20][20];

enum
	{
	BIT_MM = 0x00,
	BIT_DM = 0x01,
	BIT_IM = 0x02,
	BIT_xM = 0x03,

	BIT_DD = 0x00,
	BIT_MD = 0x04,
	//  ID not allowed
	BIT_xD = 0x04,

	BIT_II = 0x00,
	BIT_MI = 0x08,
	//  DI not allowed
	BIT_xI = 0x08,
	};
