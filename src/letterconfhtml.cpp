#include "muscle.h"
#include "heatmapcolors.h"

uint GetOverlap(uint Lo1, uint Hi1, uint Lo2, uint Hi2)
	{
	uint MaxLo = max(Lo1, Lo2);
	uint MinHi = min(Hi1, Hi2);
	if (MaxLo > MinHi)
		return 0;
	return MinHi - MaxLo + 1;
	}

static void GetConfRanges_Gapped(const char *Seq, uint L, 
  vector<char> &Confs, vector<uint> &Los, vector<uint> &His)
	{
	Confs.clear();
	Los.clear();
	His.clear();

	char CurrConf = Seq[0];
	uint Lo = 0;
	for (uint i = 1; i <= L; ++i)
		{
		char Conf = (i == L ? 0 : Seq[i]);
		if (Conf != CurrConf)
			{
			uint Hi = i - 1;

			uint n = SIZE(His);
			if (n > 0)
				asserta(Lo == His[n-1] + 1);

			Confs.push_back(CurrConf);
			Los.push_back(Lo);
			His.push_back(Hi);

			Lo = i;
			CurrConf = Conf;
			}
		}
	}

static void HTML_Head(FILE *f)
	{
	if (f == 0)
		return;
	fprintf(f,
"<!DOCTYPE html>\n"
"<html lang=\"en\">\n"
"<head>\n"
"    <meta charset=\"utf-8\">\n"
"    <title>Muscle5 alignment</title>\n"
"    <style>\n"
"        .MonoBold {font-family: monospace; font-weight: bold; font-size: 16px;}\n"
	);

	for (uint i = 0; i < 10; ++i)
		fprintf(f,
"        .Style%u {background-color: #%s;}\n", i, g_HeatmapColors_HTML[i]);

	fprintf(f,
"        .StyleG {background-color: #e6e6e6;}\n"
"        .StyleL {font-weight: normal; }\n"
"        }\n"
"    </style>\n"
"</head>\n"
"<body>\n"
"<span class=\"MonoBold\">\n"
"        <br />\n"
"        <br />\n"
	);
	}

static void HTML_Foot(FILE *f)
	{
	if (f == 0)
		return;
	fprintf(f,
"        <br />\n"
"        <span style=\"font-family:serif; font-weight:normal\">Confidence high</span>\n"
	);

	fprintf(f, "        ");
	for (int i = 9; i >= 0; --i)
		fprintf(f, "<span class=\"Style%c\">%c</span>", '0' + i, '0' + i);
	fprintf(f,
"        <span style=\"font-family:serif; font-weight:normal\">low</span>\n"
	);

	fprintf(f, "\n");
	fprintf(f,

//"        <span class=\"Style0\">0</span><span class=\"Style1\">1</span><span class=\"Style2\">2</span><span class=\"Style3\">3</span><span class=\"Style4\">4</span><span class=\"Style5\">5</span><span class=\"Style6\">6</span><span class=\"Style7\">7</span><span class=\"Style8\">8</span><span class=\"Style9\">9</span>\n"
"        </span>\n"
"    </body>\n"
"</html>\n"
	);
	}

void WriteLetterConfHTML(const string &FileName, const MSA &Ref, const MSA &ConfAln)
	{
	if (FileName.empty())
		return;

	FILE *fOut = CreateStdioFile(FileName);
	HTML_Head(fOut);
	const uint SeqCount = ConfAln.GetSeqCount();
	const uint ColCount = ConfAln.GetColCount();
	if (Ref.GetSeqCount() != SeqCount || Ref.GetColCount() != ColCount)
		Die("-ref has different number of rows or columns");

	vector<string> Labels;
	uint MaxLabelLength = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Label = (string) ConfAln.GetSeqName(SeqIndex);
		const string RefLabel = (string) Ref.GetSeqName(SeqIndex);
		if (Label != RefLabel)
			Die("-ref labels do not match, seq %u input=%s ref=%s",
			  SeqIndex + 1, Label.c_str(), RefLabel.c_str());
		MaxLabelLength = max(MaxLabelLength, SIZE(Label));
		Labels.push_back(Label);
		}

	const uint ROWLEN = 80;
	unsigned BlockCount = (ColCount + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned BlockLo = BlockIndex*ROWLEN;
		unsigned BlockHi = BlockLo + ROWLEN - 1;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const string &Label = Labels[SeqIndex];
			const uint n = SIZE(Label);
			fputs("<span class=\"StyleL\">", fOut);
			for (uint k = n; k < MaxLabelLength; ++k)
				fputs("&nbsp;", fOut);
			fputs(Label.c_str(), fOut);
			fputs("</span>&nbsp;&nbsp;", fOut);
			const char *Seq = ConfAln.m_szSeqs[SeqIndex];

			vector<char> Confs;
			vector<uint> Los;
			vector<uint> His;
			GetConfRanges_Gapped(Seq, ColCount, Confs, Los, His);

			const uint RangeCount = SIZE(Confs);
			asserta(SIZE(His) == RangeCount);
			asserta(SIZE(Los) == RangeCount);
			for (uint i = 0; i < RangeCount; ++i)
				{
				char Conf = Confs[i];
				uint Lo = Los[i];
				uint Hi = His[i];

				uint Overlap = GetOverlap(BlockLo, BlockHi, Lo, Hi);
				if (Overlap == 0)
					continue;

				if (i > 0)
					asserta(Lo == His[i-1] + 1);
				if (Conf >= '0' && Conf <= '9')
					fprintf(fOut, "<span class=\"Style%c\">", Conf);
				else if (Conf == '-')
					fprintf(fOut, "<span class=\"StyleG\">");
				else
					Die("Bad conf=%c", Conf);

				for (uint Pos = Lo; Pos <= Hi; ++Pos)
					{
					if (Pos < BlockLo || Pos > BlockHi)
						continue;
					char c = Ref.GetChar(SeqIndex, Pos);
					fputc(c, fOut);
					}

				fprintf(fOut, "</span>");
				}
			fprintf(fOut, "<br />\n");
			}
		fprintf(fOut, "<br />\n");
		fprintf(fOut, "<br />\n");
		}
	HTML_Foot(fOut);
	}

void cmd_letterconf_html()
	{
	extern bool g_FASTA_AllowDigits;
	g_FASTA_AllowDigits = true;
	MSA ConfAln;
	ConfAln.FromFASTAFile(opt(letterconf_html));
	g_FASTA_AllowDigits = false;

	const string &RefFileName = opt(output);
	if (RefFileName.empty())
		Die("Must set -ref");
	MSA Ref;
	Ref.FromFASTAFile_PreserveCase(opt(ref));

	const string &OutputFileName = opt(output);
	if (OutputFileName.empty())
		Die("Must set -output");

	WriteLetterConfHTML(OutputFileName, Ref, ConfAln);
	}
