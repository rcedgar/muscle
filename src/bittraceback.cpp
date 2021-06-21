#include "muscle.h"
#include "pwpath.h"

#define TRACE 0

static char XlatEdgeType(char c)
	{
	if ('E' == c)
		return 'D';
	if ('J' == c)
		return 'I';
	return c;
	}

static const char *BitsToStr(char Bits)
	{
	static char Str[] = "xM xD xI";

	switch (Bits & BIT_xM)
		{
	case BIT_MM:
		Str[0] = 'M';
		break;
	case BIT_DM:
		Str[0] = 'D';
		break;
	case BIT_IM:
		Str[0] = 'I';
		break;
		}

	switch (Bits & BIT_xD)
		{
	case BIT_MD:
		Str[3] = 'M';
		break;
	case BIT_DD:
		Str[3] = 'D';
		break;
		}

	switch (Bits & BIT_xI)
		{
	case BIT_MI:
		Str[6] = 'M';
		break;
	case BIT_II:
		Str[6] = 'I';
		break;
		}

	return Str;
	}

static inline char XChar(char Bits, char cType)
	{
	switch (cType)
		{
	case 'M':
		{
		switch (Bits & BIT_xM)
			{
		case BIT_MM:
			return 'M';
		case BIT_DM:
			return 'D';
		case BIT_IM:
			return 'I';
#if	DOUBLE_AFFINE
		case BIT_EM:
			return 'E';
		case BIT_JM:
			return 'J';
#endif
			}
		Quit("Huh!?");
		return '?';
		}
	case 'D':
		{
		switch (Bits & BIT_xD)
			{
		case BIT_MD:
			return 'M';
		case BIT_DD:
			return 'D';
			}
		Quit("Huh!?");
		return '?';
		}
	case 'I':
		{
		switch (Bits & BIT_xI)
			{
		case BIT_MI:
			return 'M';
		case BIT_II:
			return 'I';
			}
		Quit("Huh!?");
		return '?';
		}
#if	DOUBLE_AFFINE
	case 'E':
		{
		switch (Bits & BIT_xE)
			{
		case BIT_ME:
			return 'M';
		case BIT_EE:
			return 'E';
			}
		Quit("Huh!?");
		return '?';
		}
	case 'J':
		{
		switch (Bits & BIT_xJ)
			{
		case BIT_MJ:
			return 'M';
		case BIT_JJ:
			return 'J';
			}
		Quit("Huh!?");
		return '?';
		}
#endif
	default:
		Quit("Huh?");
		return '?';
		}
	}

void BitTraceBack(char **TraceBack, unsigned uLengthA, unsigned uLengthB,
  char LastEdge, PWPath &Path)
	{
#if	TRACE
	Log("BitTraceBack\n");
#endif
	Path.Clear();

	PWEdge Edge;
	Edge.uPrefixLengthA = uLengthA;
	Edge.uPrefixLengthB = uLengthB;
	char Bits = TraceBack[uLengthA][uLengthB];
	Edge.cType = LastEdge;
	for (;;)
		{
#if	TRACE
		Log("Prepend %c%d.%d\n", Edge.cType, Edge.uPrefixLengthA, Edge.uPrefixLengthB);
#endif
		char cSave = Edge.cType;
		Edge.cType = XlatEdgeType(cSave);
		Path.PrependEdge(Edge);
		Edge.cType = cSave;

		unsigned PLA = Edge.uPrefixLengthA;
		unsigned PLB = Edge.uPrefixLengthB;
		char Bits = TraceBack[PLA][PLB];
		char NextEdgeType = XChar(Bits, Edge.cType);
#if	TRACE
		Log("XChar(%s, %c) = %c\n", BitsToStr(Bits), Edge.cType, NextEdgeType);
#endif
		switch (Edge.cType)
			{
		case 'M':
			{
			if (Edge.uPrefixLengthA == 0)
				Quit("BitTraceBack MA=0");
			if (Edge.uPrefixLengthB == 0)
				Quit("BitTraceBack MA=0");
			--(Edge.uPrefixLengthA);
			--(Edge.uPrefixLengthB);
			break;
			}
		case 'D':
		case 'E':
			{
			if (Edge.uPrefixLengthA == 0)
				Quit("BitTraceBack DA=0");
			--(Edge.uPrefixLengthA);
			break;
			}
		case 'I':
		case 'J':
			{
			if (Edge.uPrefixLengthB == 0)
				Quit("BitTraceBack IB=0");
			--(Edge.uPrefixLengthB);
			break;
			}
		default:
			Quit("BitTraceBack: Invalid edge %c", Edge);
			}

		if (0 == Edge.uPrefixLengthA && 0 == Edge.uPrefixLengthB)
			break;

		Edge.cType = NextEdgeType;
		}

#if	TRACE
	Path.LogMe();
#endif
	}
