#include "muscle.h"
#include "masm.h"

void GetMegaProfileAASeq(const vector<vector<byte> > &Profile, string &Seq)
	{
	Seq.clear();
	uint PI = UINT_MAX;
	for (uint i = 0; i < SIZE(Mega::m_FeatureNames); ++i)
		{
		if (Mega::m_FeatureNames[i] == "AA")
			{
			PI = i;
			break;
			}
		}
	asserta(PI != UINT_MAX);
	const uint L = SIZE(Profile);
	for (uint i = 0; i < L; ++i)
		Seq += g_LetterToCharAmino[Profile[i][PI]];
	}

void WriteLocalAln_MASM(FILE *f, const string &LabelA, const MASM &MA,
  const string &LabelB, const vector<vector<byte> > &PB,
  uint Loi, uint Loj, const char *Path)
	{
	if (f == 0)
		return;
	string strA;
	string strB;
	MA.GetConsensusSeq(strA);
	GetMegaProfileAASeq(PB, strB);
	const byte *A = (const byte *) strA.c_str();
	const byte *B = (const byte *) strB.c_str();
	const unsigned BLOCK_SIZE = 80;
	uint ColLo = 0;
	uint ColHi = (unsigned) strlen(Path) - 1;

	asserta(ColHi >= ColLo);

	unsigned PosA = Loi;
	unsigned PosB = Loj;
	unsigned ColFrom = ColLo;
	for (;;)
		{
		if (ColFrom > ColHi)
			break;
		unsigned ColTo = ColFrom + BLOCK_SIZE - 1;
		if (ColTo > ColHi)
			ColTo = ColHi;

		fprintf(f, "\n");
		unsigned i0 = PosA;
		unsigned j0 = PosB;
		WriteARow(f, A, Path, PosA, ColFrom, ColTo, LabelA);
		WriteAnnotRow(f, A, B, Path, i0, j0, ColFrom, ColTo);
		WriteBRow(f, B, Path, PosB, ColFrom, ColTo, LabelB);

		ColFrom += BLOCK_SIZE;
		}
	}
