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
	FILE *fOut = CreateStdioFile(opt(output));

	MSA Aln;
	TextFile TF(MSAFileName.c_str());
	Aln.FromFASTAFile(TF);
	TF.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressLog("%u seqs, %u cols\n", SeqCount, ColCount);

	uint MinBlockLength = 8;
	uint MinBlockSeqs = 8;
	if (optset_min_core_block_cols)
		MinBlockLength = opt(min_core_block_cols);
	if (optset_min_core_block_seqs)
		MinBlockSeqs = opt(min_core_block_seqs);

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
	fprintf(fOut, "core_blocks\t%u\n", BlockCount);
	for (uint BlockIdx = 0; BlockIdx < BlockCount; ++BlockIdx)
		{
		const Rect &r = Blocks[BlockIdx];
		fprintf(fOut, "block\t%u\t%u\t%u\n", BlockIdx, r.width, r.height);
		for (int SeqIdx = r.top; SeqIdx < r.top + r.height; ++SeqIdx)
			{
			string Label;
			Aln.GetSeqLabel(SeqIdx, Label);
			vector<uint> ColToPos;
			Aln.GetColToPos(SeqIdx, ColToPos);
			for (int ColIdx = r.left; ColIdx < r.left + r.width; ++ColIdx)
				fprintf(fOut, "%c", Aln.GetChar(SeqIdx, ColIdx));
			fprintf(fOut, "\t%d", ColToPos[r.left]);
			fprintf(fOut, "\t%s", Label.c_str());
			fprintf(fOut, "\n");
			}
		}
	CloseStdioFile(fOut);
	}
