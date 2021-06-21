// enums.h
// Define enum types.
// Exploit macro hacks to avoid lots of repetetive typing.
// Generally I am opposed to macro hacks because of the
// highly obscure code that results, but in this case it
// makes maintenance much easier and less error-prone.
// The idea is that this file can be included in different
// places with different definitions of s (Start), c (Case)
// and e (End). See types.h.

s(ALPHA)
c(ALPHA, Amino)
c(ALPHA, DNA)
c(ALPHA, RNA)
e(ALPHA)

s(SEQTYPE)
c(SEQTYPE, Protein)
c(SEQTYPE, DNA)
c(SEQTYPE, RNA)
c(SEQTYPE, Auto)
e(SEQTYPE)

s(ROOT)
c(ROOT, FromClustering)
c(ROOT, MidLongestSpan)
c(ROOT, MinAvgLeafDist)
e(ROOT)

s(CLUSTER)
c(CLUSTER, UPGMA)
c(CLUSTER, UPGMAMax)
c(CLUSTER, UPGMAMin)
c(CLUSTER, UPGMB)
c(CLUSTER, NeighborJoining)
e(CLUSTER)

s(JOIN)
c(JOIN, NearestNeighbor)
c(JOIN, NeighborJoining)
e(JOIN)

s(LINKAGE)
c(LINKAGE, Min)
c(LINKAGE, Avg)
c(LINKAGE, Max)
c(LINKAGE, NeighborJoining)
c(LINKAGE, Biased)
e(LINKAGE)

s(DISTANCE)
c(DISTANCE, Kmer6_6)
c(DISTANCE, Kmer20_3)
c(DISTANCE, Kmer20_4)
c(DISTANCE, Kbit20_3)
c(DISTANCE, Kmer4_6)
c(DISTANCE, PctIdKimura)
c(DISTANCE, PctIdLog)
c(DISTANCE, PWKimura)
c(DISTANCE, PWScoreDist)
c(DISTANCE, ScoreDist)
c(DISTANCE, MyDist)
c(DISTANCE, Edit)
e(DISTANCE)

s(PPSCORE)
c(PPSCORE, LE)
c(PPSCORE, SP)
c(PPSCORE, SV)
c(PPSCORE, SPN)
e(PPSCORE)

s(SEQWEIGHT)
c(SEQWEIGHT, None)
c(SEQWEIGHT, Henikoff)
c(SEQWEIGHT, HenikoffPB)
c(SEQWEIGHT, GSC)
c(SEQWEIGHT, ClustalW)
c(SEQWEIGHT, ThreeWay)
e(SEQWEIGHT)

s(OBJSCORE)
c(OBJSCORE, SP)				// Sum of Pairs of sequences
c(OBJSCORE, DP)				// Dynamic Programming score
c(OBJSCORE, XP)				// Cross Pairs = sum of pairs between two MSAs
c(OBJSCORE, PS)				// sum of Prof-Seq score for all seqs in MSA
c(OBJSCORE, SPF)			// sum of pairs, fast approximation
c(OBJSCORE, SPM)			// sp if <= 100 seqs, spf otherwise
e(OBJSCORE)

s(TERMGAPS)
c(TERMGAPS, Full)
c(TERMGAPS, Half)
c(TERMGAPS, Ext)
e(TERMGAPS)

#undef s
#undef c
#undef e
