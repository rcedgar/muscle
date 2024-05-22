#include "muscle.h"

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
			}
		asserta(false);
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
		asserta(false);
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
		asserta(false);
		return '?';
		}
		}
	asserta(false);
	return '?';
	}

void BitTraceBack(char **TraceBack, uint uLengthA, uint uLengthB,
  char LastEdge, string &Path)
	{
	Path.resize(0);

	uint PLA = uLengthA;
	uint PLB = uLengthB;
	char Bits = TraceBack[uLengthA][uLengthB];
	char cType = LastEdge;
	for (;;)
		{
		Path.push_back(cType);

		char Bits = TraceBack[PLA][PLB];
		char NextEdgeType = XChar(Bits, cType);
		switch (cType)
			{
		case 'M':
			{
			if (PLA == 0)
				Die("BitTraceBack MA=0");
			if (PLB == 0)
				Die("BitTraceBack MA=0");
			--PLA;
			--PLB;
			break;
			}
		case 'D':
			{
			if (PLA == 0)
				Die("BitTraceBack DA=0");
			--PLA;
			break;
			}
		case 'I':
			{
			if (PLB == 0)
				Die("BitTraceBack IB=0");
			--PLB;
			break;
			}
		default:
			Die("BitTraceBack: Invalid edge %c", cType);
			}

		if (0 == PLA && 0 == PLB)
			break;

		cType = NextEdgeType;
		}
	reverse(Path.begin(), Path.end());
	}
