#include "muscle.h"
#include "myutils.h"

void SetProbconsParams();
void ProcessMuscleOptions();

int main(int argc, char **argv)
	{
	MyCmdLine(argc, argv);
	if (!opt(quiet))
		{
		PrintBanner(stderr);
		if (argc < 2)
			return 0;
		}

	if (optset_output)
		{
		opt_out = opt(output);
		optset_out = true;
		}
	if (optset_out)
		{
		opt_output = opt(out);
		optset_output = true;
		}

	SetLogFileName(opt(log));
	LogProgramInfoAndCmdLine();
	SetStartTime();

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

	SetNewHandler();
	ProcessMuscleOptions();
	SetParams();
	DoMuscle();

	exit(EXIT_Success);
	}
