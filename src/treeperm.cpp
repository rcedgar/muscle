#include "muscle.h"
#include "treeperm.h"

TREEPERM StrToTREEPERM(const string &s)
	{
	if (s == "none") return TP_None;
	if (s == "abc") return TP_ABC;
	if (s == "acb") return TP_ACB;
	if (s == "bca") return TP_BCA;
	if (s == "all") return TP_All;
	Die("Invalid perm '%s'", s.c_str());
	return TP_None;
	}

const char *TREEPERMToStr(TREEPERM TP)
	{
	switch (TP)
		{
	case TP_None:	return "none";
	case TP_ABC:	return "abc";
	case TP_ACB:	return "acb";
	case TP_BCA:	return "bca";
	case TP_All:	return "all";
	default:
		break;
		}
	asserta(false);
	return "?";
	}
