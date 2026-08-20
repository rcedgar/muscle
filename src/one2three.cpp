#include "myutils.h"
#include <map>

static map<char, string> g_aa2AAA;
static map<string, char> g_AAA2aa;

static void One2Three(char aa, const string &AAA)
	{
	g_aa2AAA[aa] = AAA;
	g_AAA2aa[AAA] = aa;
	}

static bool InitOne2Three()
	{
	One2Three('A', "ALA");
	One2Three('R', "ARG");
	One2Three('N', "ASN");
	One2Three('D', "ASP");
	One2Three('B', "ASX");
	One2Three('C', "CYS");
	One2Three('Q', "GLN");
	One2Three('E', "GLU");
	One2Three('Z', "GLX");
	One2Three('G', "GLY");
	One2Three('H', "HIS");
	One2Three('I', "ILE");
	One2Three('L', "LEU");
	One2Three('K', "LYS");
	One2Three('M', "MET");
	One2Three('F', "PHE");
	One2Three('P', "PRO");
	One2Three('S', "SER");
	One2Three('T', "THR");
	One2Three('W', "TRP");
	One2Three('Y', "TYR");
	One2Three('X', "UNK");
	One2Three('V', "VAL");
	return true;
	}

static bool g_InitOne2ThreeDone = InitOne2Three();

char GetOneFromThree(const string &AAA)
	{
	map<string, char>::const_iterator p = g_AAA2aa.find(AAA);
	if (p == g_AAA2aa.end())
		return 'X';
	return p->second;
	}

void GetThreeFromOne(char aa, string &AAA)
	{
	map<char, string>::const_iterator p = g_aa2AAA.find(aa);
	if (p == g_aa2AAA.end())
		AAA = "UNK";
	else
		AAA = p->second;
	}
