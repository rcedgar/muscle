#pragma once

enum TREEPERM
	{
	TP_None = 0,
	TP_ABC = 1,
	TP_ACB = 2,
	TP_BCA = 3,
	TP_All = 4
	};

TREEPERM StrToTREEPERM(const string &s);
const char *TREEPERMToStr(TREEPERM TP);
