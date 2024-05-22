#include "muscle.h"
#include "msa.h"
#include <math.h>
#include "xdpmem.h"
#include "pathinfo.h"
#include "objmgr.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);

// Kimura distance, calculated as:
//
//		log_e(1 - p - p*p/5)
//
// where p is the fraction of residues that differ, i.e.:
//
//		p = (1 - fractional_conservation)
//
// This measure is infinite for p = 0.8541 and is considered
// unreliable for p>0.75, i.e. pct identity <25% (twilight zone).

// ClustalW uses a table lookup for values > 0.75.
// The following table was copied from the ClustalW file dayhoff.h.

static int dayhoff_pams[]={
  195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
  196,   /* 75.1% observed d; 196 PAMs estimated */
                  197,    198,    199,    200,    200,    201,    202,  203,    
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,    
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,    
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,    
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,    
  246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */          
                                  252,    253,    254,    255,    257,  258,    
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,    
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,    
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,    
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,    
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,    
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,    
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,    
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,    
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,    
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,    
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,    
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
         /* 92.9% observed; 945 PAMs */    
  988    /* 93.0% observed; 988 PAMs */
};
static int iTableEntries = sizeof(dayhoff_pams)/sizeof(dayhoff_pams[0]);

float GetKimuraDist(float FractId)
	{
	float p = 1 - FractId;

// Typical case: use Kimura's empirical formula
	if (p < 0.75)
		return -logf(1 - p - (p*p)/5);

// Per ClustalW, return 10.0 for anything over 93%
	if (p > 0.93)
		return 10.0f;

// If p >= 0.75, use table lookup
	assert(p <= 1 && p >= 0.75);

	int iTableIndex = (int) ((p - 0.75)*1000 + 0.5);
	if (iTableIndex < 0 || iTableIndex >= iTableEntries)
		Die("Internal error in MSADistKimura::ComputeDist");

	return dayhoff_pams[iTableIndex] / 100.0f;
	}

float GetFractId(const Sequence &Seqi, const Sequence &Seqj)
	{
	const uint ColCount = Seqi.GetLength();
	asserta(Seqj.GetLength() == ColCount);
	uint n = 0;
	uint N = 0;
	const char *si = Seqi.GetCharPtr();
	const char *sj = Seqj.GetCharPtr();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char ci = si[Col];
		char cj = sj[Col];
		if (isgap(ci) && isgap(cj))
			continue;
		++N;
		if (toupper(ci) == toupper(cj))
			++n;
		}
	if (N == 0)
		return 0;
	return float(n)/float(N);
	}

float GetFractId_Path(const Sequence &Seqi, const Sequence &Seqj,
  PathInfo &PI)
	{
	const uint ColCount = PI.GetColCount();
	const char *Path = PI.GetPath();
	const char *si = Seqi.GetCharPtr();
	const char *sj = Seqj.GetCharPtr();
	uint Li = Seqi.GetLength();
	uint Lj = Seqj.GetLength();
	uint Posi = 0;
	uint Posj = 0;
	uint n = 0;
	uint N = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		switch (Path[Col])
			{
		case 'M':
			{
			char ci = si[Posi++];
			char cj = sj[Posj++];
			if (isgap(ci) && isgap(cj))
				continue;
			++N;
			if (toupper(ci) == toupper(cj))
				++n;
			break;
			}

		case 'D': ++Posi; continue;
		case 'I': ++Posj; continue;
		default: asserta(false);
			}
		}
	asserta(Posi == Li);
	asserta(Posj == Lj);
	if (N == 0)
		return 0;
	return float(n)/float(N);
	}

void GetKimuraDistMx(const MultiSequence &MSA,
  vector<vector<float> > &DistMx)
	{
	const uint SeqCount = MSA.GetSeqCount();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		DistMx[i].resize(SeqCount, FLT_MAX);

	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i][i] = 0;
		const Sequence &Seqi = *MSA.GetSequence(i);
		for (uint j = 0; j < i; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Kimura dist mx");
			const Sequence &Seqj = *MSA.GetSequence(j);
			float FractId = GetFractId(Seqi, Seqj);
			float d = GetKimuraDist(FractId);
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}

void GetKimuraDistMx_Viterbi(const MultiSequence &MS,
  vector<vector<float> > &DistMx)
	{
	const uint SeqCount = MS.GetSeqCount();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		DistMx[i].resize(SeqCount, FLT_MAX);

	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	XDPMem Mem;
	PathInfo *PI = ObjMgr::GetPathInfo();
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i][i] = 0;
		const Sequence &Seqi = *MS.GetSequence(i);
		const byte *ByteSeqi = Seqi.GetBytePtr();
		uint Li = Seqi.GetLength();
		for (uint j = 0; j < i; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Kimura dist mx (pair-wise aln)");
			const Sequence &Seqj = *MS.GetSequence(j);
			const byte *ByteSeqj = Seqj.GetBytePtr();
			uint Lj = Seqj.GetLength();
			ViterbiFastMem(Mem, ByteSeqi, Li, ByteSeqj, Lj, *PI);
			float FractId = GetFractId_Path(Seqi, Seqj, *PI);
			//LogAln(ByteSeqi, Li, ByteSeqj, Lj, *PI);
			//Log(" fractid=%.3g\n", FractId);
			float d = GetKimuraDist(FractId);
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}
