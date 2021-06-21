#include "muscle.h"
#include "msa.h"
#include <math.h>

// "Standard" NJ distance: the Kimura measure.
// This is defined to be:
//
//		log_e(1 - p - p*p/5)
//
// where p is the fraction of residues that differ, i.e.:
//
//		p = (1 - fractional_conservation)
//
// This measure is infinite for p = 0.8541 and is considered
// unreliable for p >= 0.75 (according to the ClustalW docs).
// ClustalW uses a table lookup for values > 0.75.
// The following table was copied from the ClustalW file dayhoff.h.

static int dayhoff_pams[]={
  195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
  196,   /* 75.1% observed d; 196 PAMs estimated */
                  197,    198,    199,    200,    200,    201,    202,  203,    
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,    
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,    
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,    
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,    
  246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */          
                                  252,    253,    254,    255,    257,  258,    
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,    
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,    
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,    
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,    
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,    
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,    
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,    
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,    
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,    
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,    
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,    
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
         /* 92.9% observed; 945 PAMs */    
  988    /* 93.0% observed; 988 PAMs */
};
static int iTableEntries = sizeof(dayhoff_pams)/sizeof(dayhoff_pams[0]);

double KimuraDist(double dPctId)
	{
	double p = 1 - dPctId;
// Typical case: use Kimura's empirical formula
	if (p < 0.75)
		return -log(1 - p - (p*p)/5);

// Per ClustalW, return 10.0 for anything over 93%
	if (p > 0.93)
		return 10.0;

// If p >= 0.75, use table lookup
	assert(p <= 1 && p >= 0.75);
// Thanks for Michael Hoel for pointing out a bug
// in the table index calculation in versions <= 3.52.
	int iTableIndex = (int) ((p - 0.75)*1000 + 0.5);
	if (iTableIndex < 0 || iTableIndex >= iTableEntries)
		Quit("Internal error in MSADistKimura::ComputeDist");

	return dayhoff_pams[iTableIndex] / 100.0;
	}

//double MSADistKimura::ComputeDist(const MSA &msa, unsigned uSeqIndex1,
//  unsigned uSeqIndex2)
//	{
//	double dPctId = msa.GetPctIdentityPair(uSeqIndex1, uSeqIndex2);
//	return KimuraDist(dPctId);
//	}

double KimuraDistToPctId(double dKimuraDist)
	{
// Solve quadratic equation
	const double a = 0.2;
	const double b = 1;
	const double c = 1.0 - exp(-dKimuraDist);
	const double p = (-b + sqrt(b*b + 4*a*c))/(2*a);
	return 1 - p;
	}

double PctIdToHeightKimura(double dPctId)
	{
	return KimuraDist(dPctId);
	}
