#ifndef peaker_h
#define peaker_h

#include "myutils.h"

// Minimal Peaker facade for shared_structs (SpecGet* only).
class Peaker
	{
public:
	static double SpecGetFloat(const string &Spec, const string &Name, double Default);
	static void SpecGetStr(const string &Spec, const string &Name,
	  string &Str, const string &Default);
	static bool SpecGetBool(const string &Spec, const string &Name, bool Default);
	static uint SpecGetInt(const string &Spec, const string &Name, uint Default);
	};

#endif // peaker_h
