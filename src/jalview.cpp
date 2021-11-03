#include "muscle.h"
#include "heatmapcolors.h"

static void GetConfRanges_Ungapped(const char *Seq, uint L, 
  vector<char> &Confs, vector<uint> &Los, vector<uint> &His)
	{
	Confs.clear();
	Los.clear();
	His.clear();

	char CurrConf = Seq[0];
	uint Lo = 0;
	uint Pos = 0;
	if (!isgap(Seq[0]))
		++Pos;
	for (uint i = 1; i <= L; ++i)
		{
		char Conf = (i == L ? 0 : Seq[i]);
		if (Conf != CurrConf)
			{
			uint Hi = Pos - 1;

			uint n = SIZE(His);
			if (n > 0)
				asserta(Lo == His[n-1] + 1);

			Confs.push_back(CurrConf);
			Los.push_back(Lo);
			His.push_back(Hi);

			Lo = Pos;
			CurrConf = Conf;
			}
		if (!isgap(Conf))
			++Pos;
		}
	}

void WriteLetterConfJalView(const string &FileName,
  const MSA &Ref, const MSA &ConfAln)
	{
	if (FileName.empty())
		return;

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

	FILE *fOut = CreateStdioFile(FileName);
	for (uint i = 0; i < 10; ++i)
		fprintf(fOut, "LC%u\t%s\n", i, g_HeatmapColors_JalView[i]);

	fprintf(fOut, "STARTGROUP	Muscle5_LetterConfs\n");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Labels[SeqIndex];
		const char *Seq = ConfAln.m_szSeqs[SeqIndex];

		vector<char> Confs;
		vector<uint> Los;
		vector<uint> His;
		GetConfRanges_Ungapped(Seq, ColCount, Confs, Los, His);

		const uint RangeCount = SIZE(Confs);
		asserta(SIZE(His) == RangeCount);
		asserta(SIZE(Los) == RangeCount);
		for (uint i = 0; i < RangeCount; ++i)
			{
			char Conf = Confs[i];
			uint Lo = Los[i];
			uint Hi = His[i];
			if (Conf == '-')
				continue;
			fprintf(fOut, "-	%s	%u	%u	%u	LC%c\n",
			  Label.c_str(), SeqIndex, Lo+1, Hi+1, Conf);
			}
		}
	fprintf(fOut, "ENDGROUP	Muscle5_LetterConfs\n");
	}
