#include "muscle.h"
#include "mega.h"
#include "mpcflat_mega.h"

void CheckMegaOpts(bool Nucleo)
	{
	if (optset_diversified)
		Die("-diversified not supported");
	if (optset_replicates)
		Die("-replicates not supported");
	if (optset_stratified)
		Die("-stratified not supported");
	if (optset_perm)
		Die("-perm not supported");
	if (optset_replicates)
		Die("-replicates not supported");
	if (opt_output.find('@') != string::npos)
		Die("-output pattern not supported");
	if (optset_nt)
		Die("proteins required, -nt not supported");
	if (Nucleo)
		Die("input is nucleotides, proteins required");
	}

static void Align_Mega(MPCFlat_mega &M, MultiSequence &InputSeqs,
  uint PerturbSeed, TREEPERM TP, FILE *fOut)
	{
	if (fOut == 0)
		return;

   bool Nucleo = (g_Alpha == ALPHA_Nucleo);
	HMMParams HP;
	if (optset_hmmin)
		HP.FromFile(opt(hmmin));
	else
		HP.FromDefaults(Nucleo);
	HP.CmdLineUpdate();
	if (PerturbSeed > 0)
		{
		ResetRand(PerturbSeed);
		HP.PerturbProbs(PerturbSeed);
		}
	HP.ToPairHMM();

	M.m_TreePerm = TP;
	M.Run(&InputSeqs);

	asserta(M.m_MSA != 0);
	M.m_MSA->WriteMFA(fOut);
	}

void cmd_align_mega()
	{
	Mega MM;
	MM.FromFile(g_Arg1);

	MultiSequence InputSeqs;
	InputSeqs.FromStrings(MM.m_Labels, MM.m_Seqs);
	const uint InputSeqCount = InputSeqs.GetSeqCount();

	ShowSeqStats(InputSeqs);
	CheckMegaOpts(false);

	FILE *fOut = CreateStdioFile(opt(output));

	SetAlpha(ALPHA_Amino);

	MPCFlat_mega M(MM);
	
	if (optset_consiters)
		M.m_ConsistencyIterCount = opt(consiters);
	if (optset_refineiters)
		M.m_RefineIterCount = opt(refineiters);

	uint PerturbSeed = 0;
	if (optset_perturb)
		PerturbSeed = opt(perturb);

	TREEPERM TP = TP_None;
	fOut = CreateStdioFile(opt(output));
	Align_Mega(M, InputSeqs, PerturbSeed, TP, fOut);
	CloseStdioFile(fOut);
	}
