#ifndef types_h
#define types_h

typedef unsigned char byte;
// typedef unsigned int ushort;

typedef float SCOREMATRIX[32][32];
typedef SCOREMATRIX *PTR_SCOREMATRIX;

class MSA;
class Seq;
class ClusterTree;
class DistFunc;
class TextFile;
class PWPath;
class Tree;
class SeqVect;
class DistCalc;

struct ProgNode;
struct ProfPos;

#if	SINGLE_AFFINE
// Compress M, D and I trace-back matrices into 4 bits
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

#endif

#if	DOUBLE_AFFINE
// Compress M, D, E, I and J trace-back matrices into 7 bits
enum
	{
	BIT_MM = 0x00,
	BIT_DM = 0x01,
	BIT_EM = 0x02,
	BIT_IM = 0x03,
	BIT_JM = 0x04,
	BIT_xM = 0x07,

	BIT_DD = 0x00,
	BIT_MD = 0x08,
	// [EIJ]D not sallowed
	BIT_xD = 0x08,

	BIT_EE = 0x00,
	BIT_ME = 0x10,
	//  [DDJ]E not allowed
	BIT_xE = 0x10,

	BIT_II = 0x00,
	BIT_MI = 0x20,
	//  [EDJ]I not allowed
	BIT_xI = 0x20,

	BIT_JJ = 0x00,
	BIT_MJ = 0x40,
	//  [EDI]J not allowed
	BIT_xJ = 0x40,
	};
#endif

enum EXIT
	{
	EXIT_Success = 0,
	EXIT_NotStarted = 1,
	EXIT_FatalError = 2,
	EXIT_Except = 3,
	};

enum NODECMP
	{
	NODECMP_Undefined = 0,
	NODECMP_Same = 0,		// equivalent to node in old tree
	NODECMP_Diff = 1,		// equivalent & parent is changed
	NODECMP_Changed = 2		// no equivalent node in old tree
	};

// Declare enums using macro hacks (see enums.h).
#define s(t)	enum t { t##_Undefined = 0,
#define c(t, x)	t##_##x,
#define e(t)	};
#include "enums.h"

// Declare conversion function XXXToStr(XXX x)
// for each enum type XXX.
#define	s(t)	const char *t##ToStr(t x);
#define c(t, x)	/* empty */
#define e(t)	/* empty */
#include "enums.h"

// Declare conversion function StrToXXX(const char *Str)
// for each enum type XXX.
#define	s(t)	t StrTo##t(const char *Str);
#define c(t, x)	/* empty */
#define e(t)	/* empty */
#include "enums.h"

const char *BoolToStr(bool b);
const char *SecsToStr(unsigned long Secs);

#endif // types_h
