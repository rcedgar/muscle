#include "myutils.h"
#include "muscle.h"

struct VALUE_OPT
	{
	const char *m_pstrName;
	const char *m_pstrValue;
	};

struct FLAG_OPT
	{
	const char *m_pstrName;
	bool m_bSet;
	};

static VALUE_OPT ValueOpts[] =
	{
	"in",				0,
	"in1",				0,
	"in2",				0,
	"out",				0,
	"MaxIters",			0,
	"MaxHours",			0,
	"GapOpen",			0,
	"GapOpen2",			0,
	"GapExtend",		0,
	"GapExtend2",		0,
	"GapAmbig",			0,
	"Center",			0,
	"SmoothScoreCeil",	0,
	"MinBestColScore",	0,
	"MinSmoothScore",	0,
	"ObjScore",			0,
	"SmoothWindow",		0,
	"RefineWindow",		0,
	"FromWindow",		0,
	"ToWindow",			0,
	"SaveWindow",		0,
	"WindowOffset",		0,
	"FirstWindow",		0,
	"AnchorSpacing",	0,
	"Log",				0,
	"LogA",				0,
	"MaxTrees",			0,
	"SUEFF",			0,
	"Distance",			0,
	"Distance1",		0,
	"Distance2",		0,
	"Weight",			0,
	"Weight1",			0,
	"Weight2",			0,
	"Cluster",			0,
	"Cluster1",			0,
	"Cluster2",			0,
	"Root1",			0,
	"Root2",			0,
	"Tree1",			0,
	"Tree2",			0,
	"UseTree",			0,
	"UseTree_NoWarn",	0,
	"DiagLength",		0,
	"DiagMargin",		0,
	"DiagBreak",		0,
	"Hydro",			0,
	"HydroFactor",		0,
	"SPScore",			0,
	"SeqType",			0,
	"MaxMB",			0,
	"ComputeWeights",	0,
	"MaxSubFam",		0,
	"ScoreFile",		0,
	"TermGaps",			0,
	"FASTAOut",			0,
	"CLWOut",			0,
	"CLWStrictOut",		0,
	"HTMLOut",			0,
	"MSFOut",			0,
	"PHYIOut",			0,
	"PHYSOut",			0,
	"Matrix",			0,
 	"DistMx1",			0,
 	"DistMx2",			0,
	"Weight",			0,
	};
static int ValueOptCount = sizeof(ValueOpts)/sizeof(ValueOpts[0]);

static FLAG_OPT FlagOpts[] =
	{
	"LE",					false,
	"SP",					false,
	"SV",					false,
	"SPN",					false,
	"Core",					false,
	"NoCore",				false,
	"Diags1",				false,
	"Diags2",				false,
	"Diags",				false,
	"Quiet",				false,
	"MSF",					false,
	"Verbose",				false,
	"Anchors",				false,
	"NoAnchors",			false,
	"Refine",				false,
	"RefineW",				false,
	"SW",					false,
	"Profile",				false,
	"PPScore",				false,
	"ClusterOnly",			false,
	"Brenner",				false,
	"Dimer",				false,
	"clw",					false,
	"clwstrict",			false,
	"HTML",					false,
	"Version",				false,
	"Stable",				false,
	"Group",				false,
	"FASTA",				false,
	"ProfDB",				false,
	"PAS",					false,
	"PHYI",					false,
	"PHYS",					false,
	"TomHydro",				false,
	"MakeTree",				false,
	};
static int FlagOptCount = sizeof(FlagOpts)/sizeof(FlagOpts[0]);

static bool TestSetFlagOpt(const char *Arg)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!stricmp(Arg, FlagOpts[i].m_pstrName))
			{
			FlagOpts[i].m_bSet = true;
			return true;
			}
	asserta(false);
	return false;
	}

static bool TestSetValueOpt(const char *Arg, const char *Value)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!stricmp(Arg, ValueOpts[i].m_pstrName))
			{
			if (0 == Value)
				{
				fprintf(stderr, "Option -%s must have value\n", Arg);
				exit(EXIT_NotStarted);
				}
			ValueOpts[i].m_pstrValue = strsave(Value);
			return true;
			}
	asserta(false);
	return false;
	}

bool FlagOpt(const char *Name)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!stricmp(Name, FlagOpts[i].m_pstrName))
			return FlagOpts[i].m_bSet;
	Quit("FlagOpt(%s) invalid", Name);
	return false;
	}

const char *ValueOpt(const char *Name)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!stricmp(Name, ValueOpts[i].m_pstrName))
			return ValueOpts[i].m_pstrValue;
	Quit("ValueOpt(%s) invalid", Name);
	return 0;
	}

void ProcessMuscleOptions()
	{
#define f(x)	{ if (optset_##x) TestSetFlagOpt(#x); }
	f(le)
	f(sp)
	f(sv)
	f(spn)
	f(core)
	f(nocore)
	f(diags1)
	f(diags2)
	f(diags)
	f(quiet)
	f(msf)
	f(verbose)
	f(anchors)
	f(noanchors)
	f(refine)
	f(refinew)
	f(sw)
	f(profile)
	f(ppscore)
	f(clusteronly)
	f(brenner)
	f(dimer)
	f(clw)
	f(clwstrict)
	f(html)
	f(version)
	f(stable)
	f(group)
	f(fasta)
	f(profdb)
	f(phyi)
	f(phys)
	f(tomhydro)
	f(maketree)
#undef f

#define v(x)	{ if (optset_##x) TestSetValueOpt(#x, opt_##x.c_str()); }
	v(in)
	v(in1)
	v(in2)
	v(out)
	v(maxiters)
	v(maxhours)
	v(gapopen)
	v(gapopen2)
	v(gapextend)
	v(gapextend2)
	v(gapambig)
	v(center)
	v(smoothscoreceil)
	v(minbestcolscore)
	v(minsmoothscore)
	v(objscore)
	v(smoothwindow)
	v(refinewindow)
	v(fromwindow)
	v(towindow)
	v(savewindow)
	v(windowoffset)
	v(firstwindow)
	v(anchorspacing)
	v(log)
	v(loga)
	v(maxtrees)
	v(sueff)
	v(distance)
	v(distance1)
	v(distance2)
	v(weight)
	v(weight1)
	v(weight2)
	v(cluster)
	v(cluster1)
	v(cluster2)
	v(root1)
	v(root2)
	v(tree1)
	v(tree2)
	v(usetree)
	v(usetree_nowarn)
	v(diaglength)
	v(diagmargin)
	v(diagbreak)
	v(hydro)
	v(hydrofactor)
	v(spscore)
	v(seqtype)
	v(maxmb)
	v(computeweights)
	v(maxsubfam)
	v(scorefile)
	v(termgaps)
	v(fastaout)
	v(clwout)
	v(clwstrictout)
	v(htmlout)
	v(msfout)
	v(phyiout)
	v(physout)
	v(matrix)
	v(distmx1)
	v(distmx2)
	v(weight)
#undef v
	}
