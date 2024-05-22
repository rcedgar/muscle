#include "muscle.h"
#include "uclustpd.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "omplock.h"

void GetPairs(uint Count1, uint Count2, uint TargetPairCount,
  vector<uint> &Indexes1, vector<uint> &Indexes2);
float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI);
void MakeAlnRows(const byte *XSeq, uint LX,
  const byte *YSeq, uint LY, const PathInfo &PI,
  string &RowX, string &RowY);
double GetProtDist(const char *Q, const char *T, uint ColCount);
XDPMem &GetDPMem();
double GetProtDistSeqPair(const byte *Seqi, uint Li,
  const byte *Seqj, uint Lj, string *Path = 0);

double GetProtDistSeqPair(const byte *Seqi, uint Li,
  const byte *Seqj, uint Lj, string *Path)
	{
	PathInfo *PI = ObjMgr::GetPathInfo();
	XDPMem &Mem = GetDPMem();
	ViterbiFastMem(Mem, Seqi, Li, Seqj, Lj, *PI);

	string RowX, RowY;
	MakeAlnRows(Seqi, Li, Seqj, Lj, *PI, RowX, RowY);

	const uint ColCount = PI->GetColCount();
	asserta(SIZE(RowX) == ColCount);
	asserta(SIZE(RowY) == ColCount);
	double dij = ::GetProtDist(RowX.c_str(), RowY.c_str(), ColCount);
	if (Path != 0)
		PI->GetPathStr(*Path);
	ObjMgr::Down(PI);
	return dij;
	}

double GetProtDistPairFromMFA(const MultiSequence &MFA, uint i, uint j)
	{
	uint Li;
	const byte *Seqi = MFA.GetByteSeq(i, Li);

	uint Lj;
	const byte *Seqj = MFA.GetByteSeq(j, Lj);
	double d = GetProtDistSeqPair(Seqi, Li, Seqj, Li);
	return d;
	}

double GetProtDistMFAPair(const MultiSequence &MFA1,
  const MultiSequence &MFA2, uint TargetPairCount)
	{
	asserta(TargetPairCount > 0);
	const uint SeqCount1 = MFA1.GetSeqCount();
	const uint SeqCount2 = MFA2.GetSeqCount();
	vector<uint> SeqIndexes1;
	vector<uint> SeqIndexes2;
	GetPairs(SeqCount1, SeqCount2, TargetPairCount,
	  SeqIndexes1, SeqIndexes2);
	const uint PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes2) == PairCount);

	double Sum = 0;
	for (uint i = 0; i < PairCount; ++i)
		{
		uint SeqIndex1 = SeqIndexes1[i];
		uint SeqIndex2 = SeqIndexes2[i];
		uint L1, L2;
		const byte *Seq1 = MFA1.GetByteSeq(SeqIndex1, L1);
		const byte *Seq2 = MFA2.GetByteSeq(SeqIndex2, L2);
		double d = GetProtDistSeqPair(Seq1, L1, Seq2, L2);
		Sum += d;
		}
	double Avg = Sum/PairCount;
	return Avg;
	}
