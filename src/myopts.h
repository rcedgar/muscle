#ifndef MY_VERSION
#define MY_VERSION	"5.2"
#endif

#define PROGRAM_NAME	"muscle"

////////////////////
// Commands
#define C(x)	STR_OPT(x)
#include "cmds.h"
////////////////////

STR_OPT(log)
STR_OPT(output)
STR_OPT(output1)
STR_OPT(output2)
STR_OPT(input2)
STR_OPT(joinprefix)
STR_OPT(joinpaths)
STR_OPT(linkage)
STR_OPT(joins)
STR_OPT(tsvout)
STR_OPT(html)
STR_OPT(jalview)
STR_OPT(guidetreein)
STR_OPT(guidetreeout)
STR_OPT(prefix)
STR_OPT(suffix)
STR_OPT(nodes)
STR_OPT(label)
STR_OPT(labels2)
STR_OPT(savedir)
STR_OPT(db)
STR_OPT(label1)
STR_OPT(label2)
STR_OPT(subtreeout)
STR_OPT(supertreeout)
STR_OPT(refmsa)
STR_OPT(ref)
STR_OPT(refdir)
STR_OPT(testdir)
STR_OPT(testdir1);
STR_OPT(testdir2);
STR_OPT(outdir)
STR_OPT(hmmin)
STR_OPT(hmmout)
STR_OPT(report)
STR_OPT(accalnout)
STR_OPT(perm)
STR_OPT(calnout)

UNS_OPT(threads)
UNS_OPT(consiters)
UNS_OPT(refineiters)
UNS_OPT(randseed)
UNS_OPT(paircount)
UNS_OPT(n)
UNS_OPT(splitcount)
UNS_OPT(maxcoarse)
UNS_OPT(perturb)
UNS_OPT(replicates)
UNS_OPT(maxcols)

FLT_OPT(min_cons_pct)
FLT_OPT(max_gap_fract)
FLT_OPT(minea)
FLT_OPT(super5_minea1)
FLT_OPT(super4_minea1)
FLT_OPT(super4_minea2)
FLT_OPT(pctid)
FLT_OPT(perturb_var)
FLT_OPT(minconf)

FLAG_OPT(quiet)
FLAG_OPT(compilerinfo)
FLAG_OPT(right)
FLAG_OPT(scaledist)
FLAG_OPT(eadist)
FLAG_OPT(force_super4)
FLAG_OPT(force_probcons)
FLAG_OPT(allpairs)
FLAG_OPT(nt)
FLAG_OPT(amino)
FLAG_OPT(accs)
FLAG_OPT(verbose)
FLAG_OPT(basename)
FLAG_OPT(intsuffix)
FLAG_OPT(stratified)
FLAG_OPT(diversified)
FLAG_OPT(randomchaintree)

#undef FLAG_OPT
#undef UNS_OPT
#undef FLT_OPT
#undef STR_OPT
