#pragma once

enum TREEPERM
	{
	TP_None,
	TP_ABC,
	TP_ACB,
	TP_BCA,
	TP_All
	};

TREEPERM StrToTREEPERM(const string &s);
const char *TREEPERMToStr(TREEPERM TP);
