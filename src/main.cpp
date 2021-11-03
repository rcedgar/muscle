#include "muscle.h"
#include "myutils.h"

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
			void Help();
			Help();
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

	uint CmdCount = 0;
#define C(x)	if (optset_##x) ++CmdCount;
#include "cmds.h"
	if (CmdCount > 1)
		Die("More than one command specified");

#define C(x)	\
	if (optset_##x) \
		{ \
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
