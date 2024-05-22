#include "muscle.h"
#include "m3alnparams.h"

void ReadSubstMx_Letter_FromFile(const string &FileName, float Mx[20][20])
	{
	ProgressLog("Loading %s\n", BaseName(FileName.c_str()));

	FILE *f = OpenStdioFile(FileName);

	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 21);
	for (unsigned i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		const string &s = Fields[i+1];
		asserta(s.size() == 1 && s[0] == c);
		}

	for (unsigned i = 0; i < 20; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 21);
		char c = g_LetterToCharAmino[i];
		const string &s = Fields[0];
		asserta(s.size() == 1 && s[0] == c);

		for (unsigned j = 0; j < 20; ++j)
			{
			float Score = (float) StrToFloat(Fields[j+1]);
			Mx[i][j] = Score;
			}
		}

	CloseStdioFile(f);
	}
