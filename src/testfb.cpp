#include "muscle.h"
#include "timing.h"
#include "mysparsemx.h"

#if 0
void CalcFwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);
void CalcBwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);
void CalcFwdSimple(const string &X, const string &Y,
  vector<vector<vector<float> > > &Fwd);
void CalcBwdSimple(const string &X, const string &Y,
  vector<vector<vector<float> > > &Bwd);
float CalcTotalProbFlat(const float *FlatFwd, const float *FlatBwd,
  uint LX, uint LY);
void CalcPostFlat(const float *FlatFwd, const float *FlatBwd,
  uint LX, uint LY, float *Post);
float CalcAlnFlat(const float *Post, uint LX, uint LY,
  float *DPRows, char *TB, string &Path);

static uint g_PathEqCount;
static uint g_PathDiffCount;

static void CalcBwdFlat(const string &X, const string &Y, float *Flat)
	{
	const byte *pX = (const byte *) X.c_str();
	const byte *pY = (const byte *) Y.c_str();
	uint LX = SIZE(X);
	uint LY = SIZE(Y);
	CalcBwdFlat(pX, LX, pY, LY, Flat);
	}

static void CalcFwdFlat(const string &X, const string &Y, float *Flat)
	{
	const byte *pX = (const byte *) X.c_str();
	const byte *pY = (const byte *) Y.c_str();
	uint LX = SIZE(X);
	uint LY = SIZE(Y);
	CalcFwdFlat(pX, LX, pY, LY, Flat);
	}

static void LogTBA(const vector<vector<float> > &A,
  const vector<vector<char> > &TB, uint LX, uint LY)
	{
	Log("\n");
	Log("       ");
	for (uint j = 0; j <= LY; ++j)
		Log("  %8u", j);
	Log("\n");

	for (uint i = 0; i <= LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j <= LY; ++j)
			Log("  %8.3g", A[i][j]);
		Log("\n");
		}

	Log("\n");
	Log("       ");
	for (uint j = 0; j <= LY; ++j)
		Log("  %2u", j);
	Log("\n");

	for (uint i = 0; i <= LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j <= LY; ++j)
			{
			char c = TB[i][j];
			if (c == 'D')
				c = 'B';
			else if (c == 'L')
				c = 'Y';
			else if (c == 'U')
				c = 'X';
			Log("  %2c", c);
			}
		Log("\n");
		}
	}


static void LogTomPost(const vector<float> &TomPosterior, uint LX, uint LY)
	{
	Log("\n");
	Log("TomPost LX=%u LY=%u\n", LX, LY);
	Log("       ");
	for (uint j = 0; j <= LY; ++j)
		Log("  %10u", j);
	Log("\n");

	uint Ix = 0;
	for (uint i = 0; i <= LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j <= LY; ++j)
			{
			float P = TomPosterior[Ix++];
			Log("  %10.3g", P);
			}
		Log("\n");
		}
	}

static void LogFlatPost(const float *MyPost, uint LX, uint LY)
	{
	Log("\n");
	Log("MyPost LX=%u LY=%u\n", LX, LY);
	Log("       ");
	for (uint j = 0; j < LY; ++j)
		Log("  %10u", j);
	Log("\n");

	uint Ix = 0;
	for (uint i = 0; i < LX; ++i)
		{
		Log("[%3u]  ", i);
		for (uint j = 0; j < LY; ++j)
			{
			float P = MyPost[Ix++];
			Log("  %10.3g", P);
			}
		Log("\n");
		}
	}

static bool PostEq(float x, float y)
	{
	if (x == y)
		return true;
	if (abs(x) < POSTERIOR_CUTOFF && abs(x) < POSTERIOR_CUTOFF)
		return true;
    double X = fabs(x);
    double Y = fabs(y);
    double Max = max(X, Y);
    double Diff = fabs(X-Y);
    bool Same = Diff < Max*0.05;
	if (Same)
		return true;
	return false;
	}

static void CmpPost(const vector<float> &TomPosterior, const float *PostFlat, 
  uint LX, uint LY)
	{
	for (uint i = 0; i < LX; ++i)
		{
		for (uint j = 0; j < LY; ++j)
			{
			float TomP = TomPosterior[(LY+1)*(i+1) + j + 1];
			float MyP = PostFlat[LY*i + j];
			if (!PostEq(TomP, MyP))
				{
				LogTomPost(TomPosterior, LX, LY);
				LogFlatPost(PostFlat, LX, LY);
				Die("CmpPost i=%u j=%u Tom=%.5g My=%.5g",
				  i, j, TomP, MyP);
				}
			}
		}
	}

void CvtFlat(const float *Flat, uint LX, uint LY,
  vector<vector<vector<float> > > &Mxs)
	{
	const uint LY1 = LY+1;

	Mxs.clear();
	Mxs.resize(HMMSTATE_COUNT);
	for (uint s = 0; s < HMMSTATE_COUNT; ++s)
		{
		Mxs[s].resize(LX+1);
		for (uint i = 0; i <= LX; ++i)
			Mxs[s][i].resize(LY+1, INVALID_LOG);
		}

	for (uint s = 0; s < HMMSTATE_COUNT; ++s)
		for (uint i = 0; i <= LX; ++i)
			for (uint j = 0; j <= LY; ++j)
				Mxs[s][i][j] = FLATMX(s, i, j);
	}

static void Test(const string &X, const string &Y, bool DoFwd, bool DoBwd)
	{
	Sequence &SeqX = *NewSequence();
	Sequence &SeqY = *NewSequence();

	const uint LX = SIZE(X);
	const uint LY = SIZE(Y);

	float *FlatFwd = myalloc(float, (LX+1)*(LY+1)*HMMSTATE_COUNT);
	float *FlatBwd = myalloc(float, (LX+1)*(LY+1)*HMMSTATE_COUNT);
	float *PostFlat = myalloc(float, LX*LY);

	vector<vector<vector<float> > > TomMxsFwd;
	vector<vector<vector<float> > > TomMxsBwd;
	vector<vector<vector<float> > > SimpleMxsFwd;
	vector<vector<vector<float> > > SimpleMxsBwd;
	vector<vector<vector<float> > > FlatMxsFwd;
	vector<vector<vector<float> > > FlatMxsBwd;

	SeqX.FromString("X", X);
	SeqY.FromString("Y", Y);

	SetAlpha(ALPHA_Amino);
	InitProbcons();

	vector<float> *Fwd = 0;
	vector<float> *Bwd = 0;
	vector<float> *TomPosterior = 0;

	if (DoFwd)
		{
		Fwd = PairHMM::ComputeForwardMatrix(&SeqX, &SeqY);
		PairHMM::ConvertFBMxs(*Fwd, LX, LY, TomMxsFwd);
		g_Toms = &TomMxsFwd;
		CalcFwdSimple(X, Y, SimpleMxsFwd);
		CalcFwdFlat(X, Y, FlatFwd);
		CvtFlat(FlatFwd, LX, LY, FlatMxsFwd);
		CmpFBMxs("Tom-Simple-Fwd", X, Y, TomMxsFwd, SimpleMxsFwd);
		CmpFBMxs("Tom-Flat-Fwd", X, Y, TomMxsFwd, FlatMxsFwd);
		}

	if (DoBwd)
		{
		Bwd = PairHMM::ComputeBackwardMatrix(&SeqX, &SeqY);
		PairHMM::ConvertFBMxs(*Bwd, LX, LY, TomMxsBwd);
		g_Toms = &TomMxsBwd;
		CalcBwdSimple(X, Y, SimpleMxsBwd);
		CalcBwdFlat(X, Y, FlatBwd);
		CvtFlat(FlatBwd, LX, LY, FlatMxsBwd);
		CmpFBMxs("Tom-Simple-Bwd", X, Y, TomMxsBwd, SimpleMxsBwd);
		CmpFBMxs("Tom-Flat-Bwd", X, Y, TomMxsBwd, FlatMxsBwd);
		}

	if (DoFwd && DoBwd)
		{
		float TomTotal = PairHMM::ComputeTotalProbability(LX, LY, *Fwd, *Bwd);
		float MyTotal = CalcTotalProbFlat(FlatFwd, FlatBwd, LX, LY);
		asserta(feq(TomTotal, MyTotal));

		TomPosterior = PairHMM::ComputePosteriorMatrix(&SeqX, &SeqY, *Fwd, *Bwd);
		CalcPostFlat(FlatFwd, FlatBwd, LX, LY, PostFlat);

		SparseMatrix SM((int) LX, (int) LY, *TomPosterior);
		SM.LogMe();

		MySparseMx MM;
		MM.FromPost(PostFlat, LX, LY);
		MM.LogMe();

		//LogTomPost(*TomPosterior, LX, LY);
		//LogFlatPost(PostFlat, LX, LY);
		CmpPost(*TomPosterior, PostFlat, LX, LY);

		pair<vector<char>*, float> Result =
		  PairHMM::ComputeAlignment(LX, LY, *TomPosterior);
		const vector<char> &TomPathVec = *Result.first;
		float TomScore = Result.second;
		string TomPath;
		for (uint Col = 0; Col < SIZE(TomPathVec); ++Col)
			TomPath += TomPathVec[Col];

		string Path;
		float *DPRows = myalloc(float, 2*(LY+1));
		char *TB = myalloc(char, (LX+1)*(LY+1));
		float MyScore = CalcAlnFlat(PostFlat, LX, LY, DPRows, TB, Path);

		Log("Scores %.3g %.3g\n", TomScore, MyScore);
		Log("\n");
		Log("  Path  %s\n", Path.c_str());
		LogAln(X, Y, Path);
		if (Path != TomPath)
			{
			Log("\n");
			Log("MyPath  %s\n", Path.c_str());
			}
		if (TomPath == Path)
			++g_PathEqCount;
		else
			++g_PathDiffCount;

		myfree(DPRows);
		myfree(TB);
		delete TomPosterior;
		TomPosterior = 0;
		}
	
	myfree(FlatBwd);
	myfree(FlatFwd);
	myfree(PostFlat);
	if (Bwd != 0)
		delete Bwd;
	if (Fwd != 0)
		delete Fwd;
	if (TomPosterior != 0)
		delete TomPosterior;
	}

static void GetRandomSeq(string &Seq, uint MinLen, uint MaxLen)
	{
	Seq.clear();
	uint L = MinLen + randu32() % (MaxLen - MinLen);
	for (uint i = 0; i < L; ++i)
		{
		uint Letter = randu32() % 20;
		char c = g_LetterToChar[Letter];
		Seq += c;
		}
	}

static void TestTiming()
	{
	vector<string> Seqs;
	vector<Sequence *> TomSeqs;
	const uint N = 50;
	const uint MAXL = 300;
	for (uint i = 0; i < N; ++i)
		{
		string Seq;
		GetRandomSeq(Seq, MAXL/2, MAXL);
		Seqs.push_back(Seq);
		Sequence* TomSeq = NewSequence();
		TomSeq->FromString("tom", Seq);
		TomSeqs.push_back(TomSeq);
		}

	double MyTicks = 0;
	{
	TICKS t1 = GetClockTicks();
	float* Flat = myalloc(float, MAXL * MAXL * HMMSTATE_COUNT);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Flat");
		const string& Seqi = Seqs[i];
		for (uint j = 0; j < N; ++j)
			{
			const string& Seqj = Seqs[j];
			CalcFwdFlat(Seqi, Seqj, Flat);
			CalcBwdFlat(Seqi, Seqj, Flat);
			}
		}
	TICKS t2 = GetClockTicks();
	MyTicks = double(t2 - t1);
	}

	double TomTicks = 0;
	{
	TICKS t3 = GetClockTicks();
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Tom");
		Sequence* Seqi = TomSeqs[i];
		for (uint j = 0; j < N; ++j)
			{
			Sequence *Seqj = TomSeqs[j];
			vector<float> *Fwd = PairHMM::ComputeForwardMatrix(Seqi, Seqj);
			vector<float> *Bwd = PairHMM::ComputeBackwardMatrix(Seqi, Seqj);
			delete Fwd;
			delete Bwd;
			}
		}
	TICKS t4 = GetClockTicks();
	TomTicks = double(t4 - t3);
	}

	ProgressLog("Me %.3g, Tom %.3g (%.1f%%)\n", MyTicks, TomTicks, 100*MyTicks/TomTicks);
	}

static void TestShort(bool DoFwd, bool DoBwd)
	{
	Test("MQTIF", "MSIF", DoFwd, DoBwd);
	Test("GATTACA", "MQTIF", DoFwd, DoBwd);
	Test("ABC", "DEF", DoFwd, DoBwd);
	Test("LQNGSEQVENCE", "QTHERSEQVENCEINSERT", DoFwd, DoBwd);
	}

static void TestLong(uint MaxL, bool DoFwd, bool DoBwd)
	{
	vector<string> Seqs;
	for (uint i = 0; i < 10; ++i)
		{
		string Seq;
		GetRandomSeq(Seq, MaxL/2, MaxL);
		Seqs.push_back(Seq);
		}

	Seqs.push_back("LSIDGKKYDTRLVATLLWFASLVLQDHVVDRYKDAADVLITETIYALLVTFSGTVVAKHGGNASGGYLTLILNCLVQLLLLIRSNIKRCGCTIGRCLVPAIIGDDGTY");
	Seqs.push_back("LEIDISKFDKSQQMIACLFEREIMKRFGFPDDLAEIWFNCRWICSFYDPVCGVSFKSDFQMKSGVASTFITNTLFLMSVIFYFWEPSPNAFGLFGGDDSLL");
	Seqs.push_back("GKFDKSQGLLALLIEIGIMRRFGAPEDLVELWYYSHMYTLLKDVKTGVSLKVIFQRKSGDAATFIGNTLFLLFVLAYYFGFNSLALALLGGDDSLL");
	Seqs.push_back("EEIDISKYDKSQGLLALMFECKLMKRFGVMWFNQHLSSHFYSQSTGVSGMTSFQRKSGDAATFAGNTFFLMAIVADSCKIEDLDICAFSGDDSVL");
	Seqs.push_back("LEIDISKYDKSQRELALEFECKLMKYFGVPSDIVELWFNAHVLTEVYDRTTKLNALIPYQRKSGDASTFIGNTLFLMAVICDLIPVSELELALFSGDDSLL");
	const uint N = SIZE(Seqs);

	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Testing long %u %u",
		  g_PathEqCount, g_PathDiffCount);
		const string &X = Seqs[i];
		for (uint j = 0; j < N; ++j)
			{
			const string &Y = Seqs[j];
			Test(X, Y, DoFwd, DoBwd);
			}
		}
	}
#endif // 0

void cmd_testfb()
	{
	//opt(testfb);
	//Test("MQTIF", "QTIF", true, true);
	//TestShort(true, true);
	//SetAlpha(ALPHA_Amino);
	//TestLong(200, true, true);
	//ProgressLog("Eq %u diff %u\n", g_PathEqCount, g_PathDiffCount);
	}
