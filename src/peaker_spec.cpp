#include "myutils.h"
#include "peaker.h"

bool g_QueryNeighborhood = true;

double Peaker::SpecGetFloat(const string &Spec, const string &Name, double Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	return StrToFloat(s);
	}

void Peaker::SpecGetStr(const string &Spec, const string &Name,
  string &Str, const string &Default)
	{
	vector<string> Fields;
	Split(Spec, Fields, ';');
	const string NameEq = Name + "=";
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, NameEq))
			{
			vector<string> Fields2;
			Split(Field, Fields2, '=');
			if (SIZE(Fields2) != 2)
				Die("expected name=value '%s'", Field.c_str());
			Str = Fields2[1];
			return;
			}
		}
	Str = Default;
	}

bool Peaker::SpecGetBool(const string &Spec, const string &Name, bool Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	if (s == "yes")
		return true;
	else if (s == "no")
		return false;
	Die("Peaker::SpecGetBool(%s)", s.c_str());
	return false;
	}

uint Peaker::SpecGetInt(const string &Spec, const string &Name, uint Default)
	{
	string s;
	SpecGetStr(Spec, Name, s, "");
	if (s == "")
		return Default;
	return StrToUint(s);
	}
