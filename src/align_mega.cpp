#include "muscle.h"
#include "mega.h"
#include "mpcflat_mega.h"


static void Align(MPCFlat_mega &M, MultiSequence &InputSeqs,
  uint PerturbSeed, TREEPERM TP, bool WriteEfaHdr, FILE *fOut)
    {
    if (fOut == 0)
        return;

   bool Nucleo = (g_Alpha == ALPHA_Nucleo);
    HMMParams HP;
    if (optset_hmmin)
        HP.FromFile(opt(hmmin));
    else
        HP.FromDefaults(Nucleo);
    if (PerturbSeed > 0)
        {
        ResetRand(PerturbSeed);
        HP.PerturbProbs(PerturbSeed);
        }
    HP.ToPairHMM();

    M.m_TreePerm = TP;
    M.Run(&InputSeqs);

    asserta(M.m_MSA != 0);

    if (WriteEfaHdr)
        {
        const char *TPStr = TREEPERMToStr(TP);
        fprintf(fOut, "<%s.%u\n", TPStr, PerturbSeed);
        }
    M.m_MSA->WriteMFA(fOut);
    }

void cmd_align_mega()
{
    Mega MM;
    MM.FromFile(g_Arg1);
    MultiSequence InputSeqs;
    InputSeqs.FromStrings(MM.m_Labels,MM.m_Seqs);
    
    //InputSeqs.LoadMFA(opt(align), true);
    const uint InputSeqCount = InputSeqs.GetSeqCount();
    if (optset_minsuper && InputSeqCount >= opt(minsuper))
        {
        InputSeqs.Clear();
        Progress("%u seqs, running Super5 algorithm\n", InputSeqCount);
        opt_super5 = opt(align);
        optset_super5 = true;
        void cmd_super5();
        cmd_super5();
        return;
        }

    const string &OutputPattern = opt(output);
    if (OutputPattern.empty())
        Die("Must set -output");

    double MeanSeqLength = InputSeqs.GetMeanSeqLength();
    uint MaxSeqLength = InputSeqs.GetMaxSeqLength();
    uint MinSeqLength = InputSeqs.GetMinSeqLength();
    ProgressLog("Input: %u seqs, avg length %.0f, max %u, min %u\n\n",
      InputSeqCount, MeanSeqLength, MaxSeqLength, MinSeqLength);
    if (MaxSeqLength > 100000)
        Die("Sequences too long, not appropriate for global alignment");
    if (MaxSeqLength > 50000)
        Warning("Long sequences, likey to crash / not globally alignable, max allowed is ~21k");
    if (InputSeqCount > 1000)
        Warning(">1k sequences, may be slow or use excessive memory, consider using -super5");

    bool OutputWildcard = OutputPattern.find('@') != string::npos;
    FILE *fOut = 0;

    bool IsNucleo = InputSeqs.GuessIsNucleo();
    if (IsNucleo)
        SetAlpha(ALPHA_Nucleo);
    else
        SetAlpha(ALPHA_Amino);

    MPCFlat_mega M(MM);
    
    if (optset_consiters)
        M.m_ConsistencyIterCount = opt(consiters);
    if (optset_refineiters)
        M.m_RefineIterCount = opt(refineiters);

    if (opt(stratified) && opt(diversified))
        Die("Cannot set both -stratified and -diversified");

    if (opt(stratified) || opt(diversified))
        {
        if (optset_perm || optset_perturb)
            Die("Cannot set -perm or -perturb with -stratified or -diversified");
        }

    uint RepCount = 1;
    if (opt(stratified))
        RepCount = 4;
    else if (opt(diversified))
        RepCount = 100;

    if (optset_replicates)
        RepCount = opt(replicates);

    if (RepCount == 1)
        {
        uint PerturbSeed = 0;
        if (optset_perturb)
            PerturbSeed = opt(perturb);

        TREEPERM TP = TP_None;
        if (optset_perm)
            TP = StrToTREEPERM(opt(perm));
        if (TP == TP_All)
            Die("-perm all not supported, use -stratified");

        string OutputFileName;
        if (OutputWildcard)
            MakeReplicateFileName(OutputPattern, TP, PerturbSeed, OutputFileName);
        else
            OutputFileName = OutputPattern;
        fOut = CreateStdioFile(OutputFileName);
        Align(M, InputSeqs, PerturbSeed, TP, false, fOut);
        CloseStdioFile(fOut);
        return;
        }

    bool Stratified = false;
    if (opt(stratified))
        {
        Stratified = true;
        RepCount *= 4;
        if (optset_perm)
            Die("Cannot set both -perm and -stratified");
        asserta(RepCount > 0);
        }

    string OutputFileName;
    if (!OutputWildcard)
        fOut = CreateStdioFile(OutputPattern);
    for (uint RepIndex = 0; RepIndex < RepCount; ++RepIndex)
        {
        uint PerturbSeed = (Stratified ? RepIndex/4 : RepIndex);
        TREEPERM TP = (optset_perm ?
          StrToTREEPERM(opt(perm)) : TREEPERM(RepIndex%4));
        ProgressLog("Replicate %u/%u, %s.%u\n",
          RepIndex+1, RepCount, TREEPERMToStr(TP), PerturbSeed);

        if (OutputWildcard)
            {
            MakeReplicateFileName(OutputPattern, TP, PerturbSeed, OutputFileName);
            fOut = CreateStdioFile(OutputFileName);
            }
        bool WriteEfaHeader = !OutputWildcard;
        Align(M, InputSeqs, PerturbSeed, TP, WriteEfaHeader, fOut);
        if (OutputWildcard)
            CloseStdioFile(fOut);
        }

    if (!OutputWildcard)
        CloseStdioFile(fOut);
}



