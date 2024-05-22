#include "muscle.h"
#include "qscorer.h"

extern bool g_FASTA_Upper;

static double HUE = 0.5;
static double SAT = 0.95;
static double h = DBL_MAX;

static void hsv_to_rgb(double h, double s, double v,
  byte &R, byte &G, byte &B)
	{
	uint h_i = uint(h*6);
	double f = h*6 - h_i;
	double p = v * (1 - s);
	double q = v * (1 - f*s);
	double t = v * (1 - (1 - f) * s);

	double r = DBL_MAX;
	double g = DBL_MAX;
	double b = DBL_MAX;
	if (h_i == 0)
		//r, g, b = v, t, p
		{
		r = v;
		g = t;
		b = p;
		}
	else if (h_i == 1)	
		//r, g, b = q, v, p
		{
		r = q;
		g = v;
		b = p;
		}
	else if (h_i == 2)
		//r, g, b = p, v, t
		{
		r = p;
		g = v;
		b = t;
		}
	else if (h_i == 3)	
		//r, g, b = p, q, v
		{
		r = p;
		g = q;
		b = v;
		}
	else if (h_i == 4)
		//r, g, b = t, p, v
		{
		r = t;
		g = p;
		b = v;
		}
	else if (h_i == 5)	
		//r, g, b = v, p, q
		{
		r = v;
		g = p;
		b = q;
		}
	else
		asserta(false);
	R = uint(r*256);
	G = uint(g*256);
	B = uint(b*256);
	assert(R < 256);
	assert(G < 256);
	assert(B < 256);
	}

static uint GetRandomColor()
	{
	if (h == DBL_MAX)
		h = double(randu32()%1000)/1000.0;

	h += 0.618033988749895;
	h = fmod(h, 1.0);
	byte R, G, B;
	hsv_to_rgb(h, HUE, SAT, R, G, B);
	uint Color = (R << 16) | (G << 8) | B;
	assert(Color <= 0xffffff);
	return Color;
	}

static uint Darken(uint Color, double Factor)
	{
	asserta(Factor >= 0 && Factor <= 1);
	byte R = (Color >> 16);
	byte G = (Color >> 8) % 0xff;
	byte B = Color % 0xff;

	R = byte(R*Factor);
	G = byte(G*Factor);
	B = byte(B*Factor);
	uint DarkColor = (R << 16) | (G << 8) | B;
	return DarkColor;
	}

static vector<uint> g_Colors;
static uint g_PrevRandomColor = 0;
static const char *GetColor(uint RefCol)
	{
	uint HexColor = UINT_MAX;
	if (RefCol >= SIZE(g_Colors))
		{
		uint n = SIZE(g_Colors);
		uint N = RefCol + 100;
		g_Colors.resize(N);
		for (uint i = n; i < N; ++i)
			{
			if (i%4 == 0)
				{
				g_PrevRandomColor = GetRandomColor();
				g_Colors[i] = g_PrevRandomColor;
				}
			uint k = i%4;
			asserta(k < 4);
			g_Colors[i] = Darken(g_PrevRandomColor, double(4 - k)/4.0);
			}
		}
	HexColor = g_Colors[RefCol];
	static char Tmp[64];
	sprintf(Tmp, "#%06x", HexColor);
	return Tmp;
	}

void cmd_cmp_msa()
	{
	const string &TestFileName = g_Arg1;
	const string &RefFileName = opt(ref);
	FILE *f = CreateStdioFile(opt(output));

	string Name;
	GetBaseName(RefFileName, Name);

	g_FASTA_Upper = false;

	MSA Test;
	MSA Ref;
	Test.FromFASTAFile(TestFileName);
	Ref.FromFASTAFile(RefFileName);

	QScorer QS;
	QS.Run(Name, Test, Ref);

	const uint TestSeqCount = Test.GetSeqCount();
	const uint TestColCount = Test.GetColCount();
	vector<vector<uint> > TestSeqIndexToRefCols(TestSeqCount);

	fprintf(f, "<html>\n");
	fprintf(f, "<body>\n");
	fprintf(f, "<span style=\"font-size:16px\">");
	fprintf(f, "<pre>");

	uint LabelCount = SIZE(QS.m_Labels);
	asserta(SIZE(QS.m_TestSeqIndexes) == LabelCount);
	asserta(SIZE(QS.m_RefSeqIndexes) == LabelCount);

	for (uint LabelIndex = 0; LabelIndex < LabelCount; ++LabelIndex)
		{
		uint TestSeqIndex = QS.m_TestSeqIndexes[LabelIndex];
		uint RefSeqIndex = QS.m_RefSeqIndexes[LabelIndex];
		if (TestSeqIndex == UINT_MAX || RefSeqIndex == UINT_MAX)
			continue;
		const vector<uint> &PosToRefCol = QS.m_PosToRefColVec[LabelIndex];

		vector<uint> &RefCols = TestSeqIndexToRefCols[TestSeqIndex];
		RefCols.resize(TestColCount, UINT_MAX);
		const char *Label = Test.GetLabel(TestSeqIndex);
		uint Pos = 0;
		for (uint TestCol = 0; TestCol < TestColCount; ++TestCol)
			{
			char tc = Test.GetChar(TestSeqIndex, TestCol);
			if (!isgap(tc))
				{
				uint RefCol = PosToRefCol[Pos];
				char rc = Ref.GetChar(RefSeqIndex, RefCol);
				bool Ok = (isalpha(rc) && toupper(rc) == toupper(tc));
				if (!Ok)
					{
					Log("\n");
					Log("Test label >%s\n", Label);
					Log(" Ref label >%s\n", Ref.GetLabel(RefSeqIndex));
					Log("     Pos %u\n", Pos);
					Log("Test col %u\n", TestCol);
					Log(" Ref col %u\n", RefCol);
					Log("Test seq %u\n", TestSeqIndex);
					Log(" Ref seq %u\n", RefSeqIndex);
					Log("tc '%c'\n", tc);
					Log("rc '%c'\n", rc);
					Die("!Ok");
					}
				asserta(RefCol != UINT_MAX);
				asserta(RefCols[TestCol] == UINT_MAX);
				if (isupper(rc))
					RefCols[TestCol] = RefCol;
				++Pos;
				}
			}
		}

	const uint ROWLEN = 100;
	unsigned BlockCount = (TestColCount + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= TestColCount)
			To = TestColCount;
		for (uint TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
			{
			if (TestSeqIndexToRefCols[TestSeqIndex].empty())
				continue;
			fprintf(f, "   ");
			const char *Label = Test.GetLabel(TestSeqIndex);
			for (unsigned TestCol = From; TestCol < To; ++TestCol)
				{
				char tc = Test.GetChar(TestSeqIndex, TestCol);
				uint RefCol = TestSeqIndexToRefCols[TestSeqIndex][TestCol];

				const char *Color = "black";
				if (RefCol == UINT_MAX)
					fprintf(f, "<span style=\"color:gray\">%c</span>", tc);
				else
					{
					if (islower(tc))
						Color = "gray";
					else
						Color = GetColor(RefCol);
					fprintf(f, "<span style=\"color:white;background-color:%s\">%c</span>", Color, tc);
					}
				}
			for (int i = To; i < ROWLEN; ++i)
				fprintf(f, " ");
			fprintf(f, "  <span style=\"color:black\">%s   </span>\n", Label);
			}
		fprintf(f, "\n\n");
		}
	fprintf(f, "</pre>");
	fprintf(f, "</span>");
	fprintf(f, "</body>");
	fprintf(f, "</html>\n");
	}
