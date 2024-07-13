#pragma once

#define TRACE	0
#define DOTONLY	0

static inline char GetTBBitM(byte **TB, uint i, uint j)
	{
	byte c = TB[i][j];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_DM)
		return 'D';
	if (c & TRACEBITS_IM)
		return 'I';
	return 'M';
	}

static inline char GetTBBitD(byte **TB, uint i, uint j)
	{
	byte c = TB[i][j+1];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_MD)
		return 'M';
	return 'D';
	}

static inline char GetTBBitI(byte **TB, uint i, uint j)
	{
	byte c = TB[i+1][j];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_MI)
		return 'M';
	return 'I';
	}

#if TRACE

static inline char GetTBBitM(byte **TB, uint i, uint j);
static inline char GetTBBitD(byte **TB, uint i, uint j);
static inline char GetTBBitI(byte **TB, uint i, uint j);

static uint g_LA;
static uint g_LB;
static vector<vector<float> > g_TraceSub;
static vector<vector<float> > g_TraceM;
static vector<vector<float> > g_TraceD;

static void INIT_TRACE(uint LA, uint LB, byte **TB)
	{
	g_LA = LA;
	g_LB = LB;
	g_TraceSub.clear();
	g_TraceM.clear();
	g_TraceD.clear();
	g_TraceM.resize(LA+1);
	g_TraceD.resize(LA+1);
	for (uint i = 0; i <= LA; ++i)
		{
		g_TraceM[i].resize(LB+1, UNINIT);
		g_TraceD[i].resize(LB+1, UNINIT);
		}

	g_TraceSub.resize(LA+1);
	for (uint i = 0; i <= LA; ++i)
		g_TraceSub[i].resize(LB+1, UNINIT);
	for (uint i = 0; i <= LA; ++i)
		{
		for (uint j = 0; j <= LB; ++j)
			TB[i][j] = TRACEBITS_UNINIT;
		}
	}

static void TRACE_M(int i, int j, float x)
	{
	asserta(i >= 0 && i <= (int) g_LA);
	asserta(j >= -1 && j <= (int) g_LB);
	if (j == -1)
		return;
	g_TraceM[i][j] = x;
	}

static void TRACE_D(int i, int j, float x)
	{
	asserta(i >= 0 && j <= (int) g_LB);
	asserta(j >= 0 && j <= (int) g_LB);
	g_TraceD[i][j] = x;
	}

static void TRACE_Sub(int i, int j, float x)
	{
	asserta(i >= 0 && j < (int) g_LB);
	asserta(j >= 0 && j < (int) g_LB);
	g_TraceSub[i][j] = x;
	}

static void LogMx(const string &Name, const vector<vector<float> > &Mx)
	{
	Log("\n");
	Log("%s:\n", Name.c_str());
	const uint RowCount = SIZE(Mx);
	const uint ColCount = SIZE(Mx[0]);
	Log("      ");
	for (uint j = 0; j < ColCount; ++j)
		Log("     [%3u]", j);
	Log("\n");

	for (uint i = 0; i < RowCount; ++i)
		{
		Log("[%3u] ", i);
		for (uint j = 0; j < ColCount; ++j)
			{
			float x = Mx[i][j];
			if (x == MINUS_INFINITY)
				Log("  %8.8s", "-inf");
			else if (x == UNINIT)
				Log("  %8.8s", ".");
			else
				Log("  %8.3g", x);
			}
		Log("\n");
		}
	}

static void DONE_TRACE(float BestScore, uint i, uint j, byte **TB)
	{
#if DOTONLY
	Log(" _");
	for (uint j = 0; j <= g_LB; ++j)
		Log("_");
	Log("_\n");
	for (uint i = 0; i <= g_LA; ++i)
		{
		Log("| ");
		for (uint j = 0; j <= g_LB; ++j)
			{
			char c = ' ';
			if (g_TraceSub[i][j] != UNINIT) c = '.';
			if (g_TraceM[i][j] != UNINIT) c = '.';
			if (g_TraceD[i][j] != UNINIT) c = '.';
			if (TB[i][j] != TRACEBITS_UNINIT) c = '.';
			Log("%c", c);
			}
		Log(" |");
		Log("\n");
		}
	Log(" _");
	for (uint j = 0; j <= g_LB; ++j)
		Log("_");
	Log("_\n");
	Log("\n");
	return;
#endif
	Log("Best score=%.3g (%u, %u)\n", BestScore, i, j);
	LogMx("Sub", g_TraceSub);
	LogMx("M", g_TraceM);
	LogMx("D", g_TraceD);
	Log("\n");
	Log("TB:\n");
	for (uint i = 0; i < g_LA; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j < g_LB; ++j)
			{
			if (TB[i][j] == TRACEBITS_UNINIT)
				Log(" -u-");
			else
				Log(" %c%c%c",
				  GetTBBitM(TB, i, j),
				  GetTBBitD(TB, i, j),
				  GetTBBitI(TB, i, j));
			}
		Log("\n");
		}
	}

#else
#define INIT_TRACE(LA, LB, TB)	/* empty */
#define TRACE_Sub(i, j, x)	/* empty */
#define TRACE_M(i, j, x)	/* empty */
#define TRACE_D(i, j, x)	/* empty */
#define DONE_TRACE(i, j, x, TB)	/* empty */
#endif
