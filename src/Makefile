CPP = ccache g++
CPPOPTS = -fopenmp -ffast-math -msse -mfpmath=sse -O3 -DNDEBUG -c

CC = ccache gcc
CCOPTS = -fopenmp -ffast-math -msse -mfpmath=sse -O3 -DNDEBUG -c

LNK = g++
LNKOPTS = -O3 -fopenmp -pthread -lpthread -static

HDRS = \
  alpha.h \
  alpha3.h \
  best3.h \
  chainer.h \
  cmds.h \
  countsort.h \
  derep.h \
  diagbox.h \
  eacluster.h \
  ensemble.h \
  filebuffer.h \
  flatmx.h \
  gobuff.h \
  heatmapcolors.h \
  help.h \
  hmmparams.h \
  hmmscores.h \
  hmmtrans.h \
  kmerscan.h \
  locallock.h \
  mpcflat.h \
  msa.h \
  multisequence.h \
  muscle.h \
  myopts.h \
  mysparsemx.h \
  myutils.h \
  pairhmm.h \
  pprog.h \
  pwpath.h \
  qscorer.h \
  qscorer3.h \
  quarts.h \
  scoretype.h \
  seq.h \
  sequence.h \
  sort.h \
  super4.h \
  super5.h \
  textfile.h \
  timing.h \
  transaln.h \
  tree.h \
  treeperm.h \
  treesplitter.h \
  types.h \
  uclust.h \
  upgma5.h \
  usage.h \
  usorter.h \

OBJS = \
  o/addconfseq.o \
  o/align.o \
  o/alignpairflat.o \
  o/allocflat.o \
  o/alnalnsflat.o \
  o/alnmsasflat.o \
  o/alnmsasflat3.o \
  o/alpha.o \
  o/alpha3.o \
  o/assertsameseqs.o \
  o/calcalnscoresparse.o \
  o/heatmapcolors.o \
  o/jalview.o \
  o/letterconf.o \
  o/diagbox.o \
  o/buildposterior3flat.o \
  o/buildpostflat.o \
  o/bwdflat3.o \
  o/calcalnflat.o \
  o/calcalnscoreflat.o \
  o/calcposteriorflat.o \
  o/conspairflat.o \
  o/consflat.o \
  o/defaulthmmparams.o \
  o/derep.o \
  o/disperse.o \
  o/dividetree.o \
  o/eacluster.o \
  o/efabestcols.o \
  o/efabestconf.o \
  o/efaexplode.o \
  o/efastats.o \
  o/ensemble.o \
  o/fa2efa.o \
  o/fwdflat3.o \
  o/getpairs.o \
  o/getpostpairsalignedflat.o \
  o/help.o \
  o/letterconfhtml.o \
  o/logdistmx.o \
  o/logmx.o \
  o/maxcc.o \
  o/mpcflat.o \
  o/eadistmx.o \
  o/eadistmxmsas.o \
  o/fasta.o \
  o/fasta2.o \
  o/mysparsemx.o \
  o/perturbhmm.o \
  o/getconsseq.o \
  o/globalinputms.o \
  o/guidetreejoinorder.o \
  o/hmmdump.o \
  o/hmmparams.o \
  o/jointrees.o \
  o/logaln.o \
  o/main.o \
  o/msastats.o \
  o/multisequence.o \
  o/permutetree.o \
  o/pprog.o \
  o/pprog2.o \
  o/pprogt.o \
  o/probcons.o \
  o/msa.o \
  o/msa2.o \
  o/myutils.o \
  o/colscoreefa.o \
  o/qscoreefa.o \
  o/randomchaintree.o \
  o/relabel.o \
  o/resample.o \
  o/stripgappyrows.o \
  o/tree.o \
  o/tree2.o \
  o/tree4.o \
  o/treefromfile.o \
  o/treetofile.o \
  o/pairhmm.o \
  o/progalnflat.o \
  o/project.o \
  o/qscore.o \
  o/qscore2.o \
  o/qscorer.o \
  o/quarts.o \
  o/refineflat.o \
  o/relaxflat.o \
  o/seb8.o \
  o/sequence.o \
  o/setprobconsparams.o \
  o/seq.o \
  o/stripgappycols.o \
  o/super4.o \
  o/super5.o \
  o/testfb.o \
  o/testlog.o \
  o/testscoretype.o \
  o/textfile.o \
  o/totalprobflat.o \
  o/tracebackflat.o \
  o/transaln.o \
  o/transq.o \
  o/treeperm.o \
  o/treesplitter.o \
  o/treesubsetnodes.o \
  o/trimtoref.o \
  o/trimtorefefa.o \
  o/uclust.o \
  o/upgma5.o \
  o/usage.o \
  o/usorter.o \

muscle5 : o/ $(OBJS)
	$(LNK) $(LNKOPTS) $(OBJS) -o o/muscle5
	strip -d o/muscle5

o/ :
	mkdir -p o/

o/addconfseq.o : addconfseq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/addconfseq.o addconfseq.cpp

o/align.o : align.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/align.o align.cpp

o/alignpairflat.o : alignpairflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alignpairflat.o alignpairflat.cpp

o/allocflat.o : allocflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/allocflat.o allocflat.cpp

o/alnalnsflat.o : alnalnsflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alnalnsflat.o alnalnsflat.cpp

o/alnmsasflat.o : alnmsasflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alnmsasflat.o alnmsasflat.cpp

o/alnmsasflat3.o : alnmsasflat3.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alnmsasflat3.o alnmsasflat3.cpp

o/alpha.o : alpha.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alpha.o alpha.cpp

o/alpha3.o : alpha3.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alpha3.o alpha3.cpp

o/assertsameseqs.o : assertsameseqs.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/assertsameseqs.o assertsameseqs.cpp

o/calcalnscoresparse.o : calcalnscoresparse.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/calcalnscoresparse.o calcalnscoresparse.cpp

o/heatmapcolors.o : heatmapcolors.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/heatmapcolors.o heatmapcolors.cpp

o/jalview.o : jalview.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/jalview.o jalview.cpp

o/letterconf.o : letterconf.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/letterconf.o letterconf.cpp

o/diagbox.o : diagbox.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/diagbox.o diagbox.cpp

o/buildposterior3flat.o : buildposterior3flat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/buildposterior3flat.o buildposterior3flat.cpp

o/buildpostflat.o : buildpostflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/buildpostflat.o buildpostflat.cpp

o/bwdflat3.o : bwdflat3.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/bwdflat3.o bwdflat3.cpp

o/calcalnflat.o : calcalnflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/calcalnflat.o calcalnflat.cpp

o/calcalnscoreflat.o : calcalnscoreflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/calcalnscoreflat.o calcalnscoreflat.cpp

o/calcposteriorflat.o : calcposteriorflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/calcposteriorflat.o calcposteriorflat.cpp

o/conspairflat.o : conspairflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/conspairflat.o conspairflat.cpp

o/consflat.o : consflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/consflat.o consflat.cpp

o/defaulthmmparams.o : defaulthmmparams.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/defaulthmmparams.o defaulthmmparams.cpp

o/derep.o : derep.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/derep.o derep.cpp

o/disperse.o : disperse.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/disperse.o disperse.cpp

o/dividetree.o : dividetree.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/dividetree.o dividetree.cpp

o/eacluster.o : eacluster.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/eacluster.o eacluster.cpp

o/efabestcols.o : efabestcols.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/efabestcols.o efabestcols.cpp

o/efabestconf.o : efabestconf.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/efabestconf.o efabestconf.cpp

o/efaexplode.o : efaexplode.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/efaexplode.o efaexplode.cpp

o/efastats.o : efastats.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/efastats.o efastats.cpp

o/ensemble.o : ensemble.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ensemble.o ensemble.cpp

o/fa2efa.o : fa2efa.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fa2efa.o fa2efa.cpp

o/fwdflat3.o : fwdflat3.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fwdflat3.o fwdflat3.cpp

o/getpairs.o : getpairs.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/getpairs.o getpairs.cpp

o/getpostpairsalignedflat.o : getpostpairsalignedflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/getpostpairsalignedflat.o getpostpairsalignedflat.cpp

o/help.o : help.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/help.o help.cpp

o/letterconfhtml.o : letterconfhtml.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/letterconfhtml.o letterconfhtml.cpp

o/logdistmx.o : logdistmx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/logdistmx.o logdistmx.cpp

o/logmx.o : logmx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/logmx.o logmx.cpp

o/maxcc.o : maxcc.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/maxcc.o maxcc.cpp

o/mpcflat.o : mpcflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/mpcflat.o mpcflat.cpp

o/eadistmx.o : eadistmx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/eadistmx.o eadistmx.cpp

o/eadistmxmsas.o : eadistmxmsas.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/eadistmxmsas.o eadistmxmsas.cpp

o/fasta.o : fasta.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fasta.o fasta.cpp

o/fasta2.o : fasta2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fasta2.o fasta2.cpp

o/mysparsemx.o : mysparsemx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/mysparsemx.o mysparsemx.cpp

o/perturbhmm.o : perturbhmm.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/perturbhmm.o perturbhmm.cpp

o/getconsseq.o : getconsseq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/getconsseq.o getconsseq.cpp

o/globalinputms.o : globalinputms.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/globalinputms.o globalinputms.cpp

o/guidetreejoinorder.o : guidetreejoinorder.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/guidetreejoinorder.o guidetreejoinorder.cpp

o/hmmdump.o : hmmdump.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/hmmdump.o hmmdump.cpp

o/hmmparams.o : hmmparams.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/hmmparams.o hmmparams.cpp

o/jointrees.o : jointrees.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/jointrees.o jointrees.cpp

o/logaln.o : logaln.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/logaln.o logaln.cpp

o/main.o : main.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/main.o main.cpp

o/msastats.o : msastats.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/msastats.o msastats.cpp

o/multisequence.o : multisequence.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/multisequence.o multisequence.cpp

o/permutetree.o : permutetree.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/permutetree.o permutetree.cpp

o/pprog.o : pprog.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pprog.o pprog.cpp

o/pprog2.o : pprog2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pprog2.o pprog2.cpp

o/pprogt.o : pprogt.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pprogt.o pprogt.cpp

o/probcons.o : probcons.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/probcons.o probcons.cpp

o/msa.o : msa.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/msa.o msa.cpp

o/msa2.o : msa2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/msa2.o msa2.cpp

o/myutils.o : myutils.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/myutils.o myutils.cpp

o/colscoreefa.o : colscoreefa.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/colscoreefa.o colscoreefa.cpp

o/qscoreefa.o : qscoreefa.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/qscoreefa.o qscoreefa.cpp

o/randomchaintree.o : randomchaintree.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/randomchaintree.o randomchaintree.cpp

o/relabel.o : relabel.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/relabel.o relabel.cpp

o/resample.o : resample.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/resample.o resample.cpp

o/stripgappyrows.o : stripgappyrows.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/stripgappyrows.o stripgappyrows.cpp

o/tree.o : tree.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tree.o tree.cpp

o/tree2.o : tree2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tree2.o tree2.cpp

o/tree4.o : tree4.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tree4.o tree4.cpp

o/treefromfile.o : treefromfile.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/treefromfile.o treefromfile.cpp

o/treetofile.o : treetofile.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/treetofile.o treetofile.cpp

o/pairhmm.o : pairhmm.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pairhmm.o pairhmm.cpp

o/progalnflat.o : progalnflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/progalnflat.o progalnflat.cpp

o/project.o : project.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/project.o project.cpp

o/qscore.o : qscore.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/qscore.o qscore.cpp

o/qscore2.o : qscore2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/qscore2.o qscore2.cpp

o/qscorer.o : qscorer.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/qscorer.o qscorer.cpp

o/quarts.o : quarts.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/quarts.o quarts.cpp

o/refineflat.o : refineflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/refineflat.o refineflat.cpp

o/relaxflat.o : relaxflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/relaxflat.o relaxflat.cpp

o/seb8.o : seb8.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seb8.o seb8.cpp

o/sequence.o : sequence.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/sequence.o sequence.cpp

o/setprobconsparams.o : setprobconsparams.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/setprobconsparams.o setprobconsparams.cpp

o/seq.o : seq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seq.o seq.cpp

o/stripgappycols.o : stripgappycols.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/stripgappycols.o stripgappycols.cpp

o/super4.o : super4.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/super4.o super4.cpp

o/super5.o : super5.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/super5.o super5.cpp

o/testfb.o : testfb.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testfb.o testfb.cpp

o/testlog.o : testlog.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testlog.o testlog.cpp

o/testscoretype.o : testscoretype.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testscoretype.o testscoretype.cpp

o/textfile.o : textfile.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/textfile.o textfile.cpp

o/totalprobflat.o : totalprobflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/totalprobflat.o totalprobflat.cpp

o/tracebackflat.o : tracebackflat.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tracebackflat.o tracebackflat.cpp

o/transaln.o : transaln.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/transaln.o transaln.cpp

o/transq.o : transq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/transq.o transq.cpp

o/treeperm.o : treeperm.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/treeperm.o treeperm.cpp

o/treesplitter.o : treesplitter.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/treesplitter.o treesplitter.cpp

o/treesubsetnodes.o : treesubsetnodes.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/treesubsetnodes.o treesubsetnodes.cpp

o/trimtoref.o : trimtoref.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/trimtoref.o trimtoref.cpp

o/trimtorefefa.o : trimtorefefa.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/trimtorefefa.o trimtorefefa.cpp

o/uclust.o : uclust.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/uclust.o uclust.cpp

o/upgma5.o : upgma5.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/upgma5.o upgma5.cpp

o/usage.o : usage.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/usage.o usage.cpp

o/usorter.o : usorter.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/usorter.o usorter.cpp
