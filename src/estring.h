#ifndef pathsum_h
#define pathsum_h

void PathToEstrings(const PWPath &Path, int **ptresA, int **ptresB);
void EstringsToPath(const int esA[], const int esB[], PWPath &Path);
void MulEstrings(const int es1[], const int es2[], int esp[]);
void EstringOp(const int es[], const Seq &sIn, Seq &sOut);
unsigned EstringOp(const int es[], const Seq &sIn, MSA &a);
void LogEstring(const int es[]);
unsigned LengthEstring(const int es[]);
int *EstringNewCopy(const int es[]);

#endif	// pathsum_h
