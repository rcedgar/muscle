#include "muscle.h"

void WriteAnnotRow(FILE *f, const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi);
void WriteBRow(FILE *f, const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi);
void WriteARow(FILE *f, const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi);

void WriteLocalAln(FILE *f, const byte *A, const byte *B,
  uint Loi, uint Loj, const char *Path)
	{
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

		unsigned i0 = PosA;
		unsigned j0 = PosB;
		WriteARow(f, A, Path, PosA, ColFrom, ColTo);
		WriteAnnotRow(f, A, B, Path, i0, j0, ColFrom, ColTo);
		WriteBRow(f, B, Path, PosB, ColFrom, ColTo);
		fprintf(f, "\n");

		ColFrom += BLOCK_SIZE;
		}
	}
