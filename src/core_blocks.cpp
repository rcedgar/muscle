#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "seq.h"
#include "rect.h"

void greedy_rects(
	const std::vector<std::vector<bool>>& mat,
	int minW,
	int minH,
	std::vector<Rect> &out);

void cmd_core_blocks()
	{
	const string &MSAFileName = g_Arg1;

	MSA Aln;
	TextFile TF(MSAFileName.c_str());
	Aln.FromFASTAFile(TF);
	TF.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressLog("%u seqs, %u cols\n", SeqCount, ColCount);

	uint MinBlockLength = 8;
	uint MinBlockSeqs = 4;

	vector<vector<bool> > UngappedMx(SeqCount);
	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		UngappedMx[SeqIdx].resize(ColCount);

	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
			UngappedMx[SeqIdx][ColIdx] = !Aln.IsGap(SeqIdx, ColIdx);
		}

	vector<Rect> Blocks;
	greedy_rects(UngappedMx, MinBlockLength, MinBlockSeqs, Blocks);
	uint BlockCount = SIZE(Blocks);

	ProgressLog("%u blocks\n", BlockCount);
	for (uint BlockIdx = 0; BlockIdx < BlockCount; ++BlockIdx)
		{
		const Rect &r = Blocks[BlockIdx];
		Log("\nBlock %u width=%d\n", BlockIdx, r.width);
		for (uint SeqIdx = r.top; SeqIdx < r.top + r.height; ++SeqIdx)
			{
			string Label;
			Aln.GetSeqLabel(SeqIdx, Label);
			vector<uint> ColToPos;
			Aln.GetColToPos(SeqIdx, ColToPos);
			for (int ColIdx = r.left; ColIdx < r.left + r.width; ++ColIdx)
				Log("%c", Aln.GetChar(SeqIdx, ColIdx));
			Log(" %d", ColToPos[r.left]);
			Log(" %s", Label.c_str());
			Log("\n");
			}
		}
	}
