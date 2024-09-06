#ifndef myutils_h
#define myutils_h

#ifdef _MSC_VER
#if !defined(MYUTILS_CPP)
#define _SCL_SECURE_NO_WARNINGS
#endif // !defined(IN_MYUTILS_CPP)
#endif

#include <cstddef>
#include "myalloc.h"

#if defined(__x86_64__) || defined(_M_X64) || defined(__arm64__)
#define	BITS			64
#else
#define	BITS			32
#endif

#include <stdio.h>
#include <sys/types.h>
#include <string>
#include <string.h>
#include <memory.h>
#include <vector>
#include <math.h>
#include <stdarg.h>
#include <cstdlib>
#include <climits>
#include <float.h>
#include <omp.h>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <inttypes.h>
#include <thread>
#include "types.h"

// isatty()
#ifdef	WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifndef _MSC_VER
#define _stricmp	strcasecmp
#endif

//#ifndef max
//#define max(x, y) ((x) >= (y) ? (x) : (y))
//#define min(x, y) ((x) <= (y) ? (x) : (y))
//#endif

using namespace std;

#ifdef _MSC_VER
#include <crtdbg.h>
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

#if defined(_DEBUG) && !defined(DEBUG)
#define DEBUG	1
#endif

#if defined(DEBUG) && !defined(_DEBUG)
#define _DEBUG	1
#endif

#ifndef NDEBUG
#define	DEBUG	1
#define	_DEBUG	1
#endif

#define byte __mybyte__

typedef unsigned char byte;
typedef unsigned short uint16;
typedef unsigned uint32;
typedef unsigned uint;

// typedefs for int64 and uint64
#if		defined(_MSC_VER)
typedef __int64 int64;
typedef unsigned __int64 uint64;
#elif defined(__GNUC__)
typedef long int64;
typedef unsigned long uint64;
#else	
#error	"int64 typedefs"
#endif

#if	BITS==32
typedef uint32 uintb;
#else
typedef uint64 uintb;
#endif

#ifndef UINT32_MAX
const uint32 UINT32_MAX = (~(uint32(0)));
#endif

#ifndef UINT64_MAX
const uint64 UINT64_MAX = (~(uint64(0)));
#endif

#ifndef SIZE_T_MAX
const size_t SIZE_T_MAX = (~(size_t(1)));
#endif

void myassertfail(const char *Exp, const char *File, unsigned Line);
#undef  assert
#ifdef  NDEBUG
#define assert(exp)     ((void)0)
#define myassert(exp)     ((void)0)
#else
#define assert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#define myassert(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )
#endif
#define asserta(exp) (void)( (exp) || (myassertfail(#exp, __FILE__, __LINE__), 0) )

#define NotUsed(v)	((void *) &v)

// pom=plus or minus, tof=true or false, yon=yes or no
static inline char pom(bool Plus)	{ return Plus ? '+' : '-'; }
static inline char tof(bool x)		{ return x ? 'T' : 'F';	}
static inline char yon(bool x)		{ return x ? 'Y' : 'N';	}
static inline const char *YesOrNo(bool x)	{ return x ? "Yes" : "No"; }
static inline const char *plurals(unsigned n) { return n == 1 ? "" : "s"; }

const char *GetPlatform();
unsigned GetElapsedSecs();
void mysleep(unsigned ms);
void mylistdir(const string &DirName, vector<string> &FileNames);
void PrintVersion(FILE *f);

#define SIZE(c)	unsigned((c).size())
#define RoundUp(Bytes, BlockSize)	((Bytes) + ((BlockSize) - (Bytes)%(BlockSize)))

bool myisatty(int fd);

#ifdef _MSC_VER
#define off_t	__int64
#endif

void GetVersionString(string &s);

// Stdio functions without "nr of bytes" arg.
FILE *OpenStdioFile(const string &FileName);
FILE *CreateStdioFile(const string &FileName);
void CloseStdioFile(FILE *f);
bool ReadLineStdioFile(FILE *f, string &Line);
void ReadTabbedLine(FILE *f, vector<string> &Fields, unsigned FieldCount);
void FlushStdioFile(FILE *f);
bool StdioFileExists(const string &FileName);
void LogStdioFileState(FILE *f);
void RenameStdioFile(const string &FileNameFrom, const string &FileNameTo);
void MoveStdioFile(const string &FileName1, const string &FileName2);
void DeleteStdioFile(const string &FileName);
void WriteStdioFileStr(FILE *f, const char *s);
void Pr(FILE *f, const char *Format, ...);
void ParseFileName(const string &FileName, string &Path, string &Name);
// void ReadDir(const string &DirName, vector<string> &FileNames);

// Stdio functions with size args:
byte *ReadAllStdioFile32(FILE *f, uint32 &FileSize);
byte *ReadAllStdioFile64(FILE *f, uint64 &FileSize);

byte *ReadAllStdioFile(FILE *f, uint32 &FileSize);
byte *ReadAllStdioFile64(FILE *f, uint64 &FileSize);

byte *ReadAllStdioFile32(const string &FileName, uint32 &FileSize);
byte *ReadAllStdioFile64(const string &FileName, uint64 &FileSize);

bool ReadLineStdioFile(FILE *f, char *Line, uint32 Bytes);
bool ReadLineStdioFile64(FILE *f, char *Line, uint64 Bytes);

void SetStdioFilePos(FILE *f, uint32 Pos);
void SetStdioFilePos64(FILE *f, uint64 Pos);

uint32 GetStdioFilePos32(FILE *f);
uint64 GetStdioFilePos64(FILE *f);

uint32 GetStdioFileSize32(FILE *f);
uint64 GetStdioFileSize64(FILE *f);

#if	BITS==32
#define uintB	uint32
#define GetStdioFilePosB	GetStdioFilePos32
#define GetStdioFileSizeB	GetStdioFileSize32
#define SetStdioFilePosB	SetStdioFilePos
#else
#define uintB	uint64
#define GetStdioFilePosB	GetStdioFilePos64
#define GetStdioFileSizeB	GetStdioFileSize64
#define SetStdioFilePosB	SetStdioFilePos64
#endif

uint32 ReadStdioFile_NoFail(FILE *f, void *Buffer, uint32 Bytes);

void ReadStdioFile(FILE *f, uint32 Pos, void *Buffer, uint32 Bytes);
void ReadStdioFile64(FILE *f, uint64 Pos, void *Buffer, uint64 Bytes);

void ReadStdioFile(FILE *f, void *Buffer, uint32 Bytes);
void ReadStdioFile64(FILE *f, void *Buffer, uint64 Bytes);

void WriteStdioFile(FILE *f, uint32 Pos, const void *Buffer, uint32 Bytes);
void WriteStdioFile64(FILE *f, uint64 Pos, const void *Buffer, uint64 Bytes);

void WriteStdioFile(FILE *f, const void *Buffer, uint32 Bytes);
void WriteStdioFile64(FILE *f, const void *Buffer, uint64 Bytes);

#define MAGIC(a, b, c, d)	uint32(uint32(a)<<24 | uint32(b)<<16 | uint32(c)<<8 | (d))

void Pf(FILE *f, const char *szFormat, ...);

void Ps(string &Str, const char *szFormat, ...);
void Psa(string &Str, const char *szFormat, ...);
void Psasc(string &Str, const char *szFormat, ...);

void SetLogFileName(const string &FileName);
void Log(const char *szFormat, ...);

#if	DEBUG
void Die_(const char *szFormat, ...);
void Warning_(const char *szFormat, ...);

typedef void (*PTR_PRINTFLIKE_FN)(const char *Format, ...);

static inline PTR_PRINTFLIKE_FN DiePtr(const char *FileName, unsigned LineNr)
	{
	fprintf(stderr, "\n\n%s(%u): ", FileName, LineNr);
	Log("\n\n%s(%u): ", FileName, LineNr);
	return Die_;
	}

static inline PTR_PRINTFLIKE_FN WarningPtr(const char *FileName, unsigned LineNr)
	{
	fprintf(stderr, "\n\n%s(%u): ", FileName, LineNr);
	Log("\n\n%s(%u): ", FileName, LineNr);
	return Warning_;
	}

#define Die		(*DiePtr(__FILE__, __LINE__))
#define Warning	(*WarningPtr(__FILE__, __LINE__))
#else
void Die_(const char *szFormat, ...);
void Warning_(const char *szFormat, ...);
#define Die	Die_
#define Warning Warning_
#endif

typedef const char *(*FN_PROGRESS_CALLBACK)();
void SetPCB(FN_PROGRESS_CALLBACK PCB);
void ProgressCallback(unsigned i, unsigned N);

bool ProgressPrefix(bool On);
void ProgressStep(unsigned i, unsigned N, const char *Format, ...);
void ProgressStep64(uint64 i, uint64 N, const char *Msg);
void Progress(const char *szFormat, ...);
void Progress(const string &Str);
void ProgressLog(const char *szFormat, ...);
void ProgressLogPrefix(const char *Format, ...);
const char *GetProgressPrefixCStr();

void ProgressFileInit(FILE *f, const char *Format = 0, ...);
void ProgressFileStep(const char *Format = 0, ...);
void ProgressFileDone(const char *Format = 0, ...);

void LogElapsedTimeAndRAM();
void LogProgramInfoAndCmdLine();
void Help();

inline unsigned ustrlen(const char *s) { return (unsigned) strlen(s); }
inline unsigned ustrlen(const string &s) { return SIZE(s); }
char *mystrsave(const char *s);
unsigned myipow(unsigned x, unsigned y);
uint64 myipow64(unsigned x, unsigned y);
static inline unsigned atou(const char *s) { return unsigned(atoi(s)); }
static inline unsigned atou(const string &s) { return unsigned(atoi(s.c_str())); }

unsigned GetThreadIndex();

double GetMemUseBytes();
double GetPeakMemUseBytes();
double GetPhysMemBytes();
double GetUsableMemBytes();

// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
inline bool feq(double x, double y)
    {
	if (x == y)
		return true;
    double X = fabs(x);
    double Y = fabs(y);
    double Max = max(X, Y);
    double Diff = fabs(X-Y);
    return Diff < Max*0.01;
    }

#define asserteq(x, y)	assert(feq(x, y))
#define assertaeq(x, y)	asserta(feq(x, y))

#define	memset_zero(a, n)	memset((a), 0, (n)*sizeof(a[0]))

void ResetRand(unsigned Seed);
unsigned randu32();
uint64 randu64();
void Split(const string &Str, vector<string> &Fields, char Sep = '\t');
void StripWhiteSpace(string &Str);
bool StartsWith(const string &s, const string &t);
bool StartsWith(const char *s, const char *t);
void ToUpper(const string &s, string &t);
void ToLower(const string &s, string &t);
void Reverse(string &s);
static inline double GetRatio(double x, double y) { if (y == 0) { asserta(x == 0); return 0; } return x/y; }
static inline double GetPct(double x, double y) { return 100.0*GetRatio(x, y); }
static inline double GetPct64(int64 x, int64 y) { return 100.0*GetRatio(double(x), double(y)); }
const char *PctToStr(double Pct);
double GetMemUseBytes();
void PrintCmdLine(FILE *f);
void PrintBanner(FILE *f);

const char *MemBytesToStr(double Bytes);
static inline const char *MemBytesToStr(uint64 Bytes) { return MemBytesToStr((double) Bytes); }
double StrToMemBytes(const string &s);
double StrToFloat(const char *s, bool StarIsDblMax = false);
double StrToFloat(const string &s, bool StarIsDblMax = false);
bool IsValidFloatStr(const char *s);
bool IsValidFloatStr(const string &s);
const char *GetElapsedTimeStr(string &s);
const char *GetMaxRAMStr(string &s);

const char *BaseName(const char *PathName);

const char *IntToStr(uint64 i);
const char *IntToStr2(uint64 i);
const char *Int64ToStr(uint64 i);
const char *FloatToStr(double d);
const char *FloatToStr(uint64 i);
const char *IntFloatToStr(double d);
const char *SecsToStr(double Secs);
bool IsUintStr(const char *s);
unsigned StrToUint(const char *s, bool StarIsUnitMax = false);
unsigned StrToUint(const string &s, bool StarIsUnitMax = false);
uint64 StrToUint64(const char *s);
uint64 StrToUint64(const string &s);
bool EndsWith(const string &s, const string &t);
bool Replace(string &s, const string &a, const string &b);

double mylog2(double x);
double mylog10(double x);

void LogInt(unsigned i, unsigned w = UINT_MAX);
void Logu(unsigned u, unsigned w, unsigned prefixspaces = 2);
void Logf(float x, unsigned w, unsigned prefixspaces = 2);
const char *SecsToHHMMSS(unsigned Secs);
unsigned GetCPUCoreCount();

void MyCmdLine(int argc, char **argv);
//void CmdLineErr(const char *Format, ...);
void GetCmdLine(string &s);

#define FLAG_OPT(Name)		extern bool opt_##Name; extern bool optset_##Name; extern bool optused_##Name;
#define UNS_OPT(Name)		extern unsigned opt_##Name; extern bool optset_##Name; extern bool optused_##Name;
#define FLT_OPT(Name)		extern double opt_##Name; extern bool optset_##Name; extern bool optused_##Name;
#define STR_OPT(Name)		extern string opt_##Name; extern bool optset_##Name; extern bool optused_##Name;
#include "myopts.h"

#define opt(Name)			(optused_##Name = true, opt_##Name)
#define optd(Name, Default)	(optset_##Name ? (optused_##Name = true, opt_##Name) : Default)

void CheckUsedOpts(bool LogAll);

extern FILE *g_fLog;

void SeqToFasta(FILE *f, const string &Seq, const string &Label);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);
void SeqToFastaRC(FILE *f, const byte *Seq, unsigned L, const char *Label);
void RevCompSeq(const byte *Seq, unsigned L, byte *RCSeq);

void LogAllocs();
unsigned GetRequestedThreadCount();
void Dirize(string &Dir);

inline char mytoupper(char c) { return c & (~0x20); }
inline char mytoupper(byte c) { return c & (~0x20); }
inline bool myislower(byte c) { return (c & 0x20) != 0; }
inline bool myislower(char c) { return (c & 0x20) != 0; }
bool AccChar(char c);
void GetAccFromLabel(const string &Label, string &Acc);
void GetBaseName(const string &PathName, string &BaseName);
void Shuffle(vector<unsigned> &v);

#define NO_TRACE		0
#define TMP_TRACE		2	// true but not 1, for grep_trace
#define REMOVEME		1

extern string g_Arg1;
typedef void fn_thread_body(uint ThreadIndex, void *ptrUserData);
void RunThreads(fn_thread_body Body, void *ptrUserData);

#ifdef _MSC_VER
#define brk(x)       if (x) __debugbreak()
#else
#define brk(x)		(0)
#endif

#endif	// myutils_h
