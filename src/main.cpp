#include "muscle.h"
#include "myutils.h"

string g_Arg1;

int main(int argc, char **argv)
	{
	for (int i = 1; i < argc; ++i)
		{
		string s = string(argv[i]);
		if (s == "-h")
			{
			void Usage(FILE *f);
			Usage(stdout);
			return 0;
			}

		if (s == "-help" || s == "--help")
			{
			void Usage(FILE *f);
			Usage(stdout);
			return 0;
			}
		}

	MyCmdLine(argc, argv);
	if (!opt(quiet))
		{
		PrintBanner(stderr);
		if (argc < 2)
			return 0;
		}

	SetLogFileName(opt(log));
	LogProgramInfoAndCmdLine();

	extern vector<string> g_Argv;
	uint n = SIZE(g_Argv);
	asserta(n > 0);
	string ShortCmdLine;
	if (n > 1)
		ShortCmdLine = g_Argv[1];
	if (n > 2)
		{
		g_Arg1 = g_Argv[2];
		ShortCmdLine += " " + g_Argv[2];
		}
	if (n > 1)
		{
		ProgressPrefix(false);
		Progress("[%s]\n", ShortCmdLine.c_str() + 1);
		ProgressPrefix(true);
		}

	uint CmdCount = 0;
#define C(x)	if (optset_##x) ++CmdCount;
#include "cmds.h"
	if (CmdCount > 1)
		Die("More than one command specified");

#define C(x)	\
	if (optset_##x) \
		{ \
		g_Arg1 = opt_##x; \
		optused_##x = true; \
		void cmd_##x(); \
		cmd_##x(); \
		CheckUsedOpts(false); \
		LogElapsedTimeAndRAM(); \
		return 0; \
		}
#include "cmds.h"
#undef C

	return 0;
	}
