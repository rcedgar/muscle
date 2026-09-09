#!/usr/bin/python3

import os
import sys
import argparse

Usage = \
(
"Convert Visual Studio .vcxproj file in current directory to Makefile or make.bash, generate git hash, and run make. Links d5 myutils."
)

D5_MYUTILS_INC = "../../d5/myutils"
D5_MYUTILS_LIB_V2 = "../../d5/build-wsl-v2/libmyutils.a"
D5_MYUTILS_LIB_V3 = "../../d5/build-wsl-v3/libmyutils.a"

AP = argparse.ArgumentParser(description = Usage)

# Value opts
AP.add_argument("--std", required=False, help="C++ standard option for GCC, e.g. c++11 or c++17 (default none)")
AP.add_argument("--deletes", required=False, help="source filenames to delete separated by +")
AP.add_argument("--cppcompiler", required=False, default="g++", help="C++ compiler command name default g++)")
AP.add_argument("--ccompiler", required=False, default="gcc", help="C compiler command name default gcc)")
AP.add_argument("--binary", help="Binary filename (path not allowed, default project name")

# Flag opts
AP.add_argument("--debug", required=False, action='store_true', help="Debug build")
AP.add_argument("--profile", required=False, action='store_true', help="Profile build")
AP.add_argument("--openmp", required=False, action='store_true', help="Requires OMP")
AP.add_argument("--santhread", required=False, action='store_true', help="Set -fsanitize=thread")
AP.add_argument("--sanaddr", required=False, action='store_true', help="Set -fsanitize=address")
AP.add_argument("--pthread", required=False, action='store_true', help="Requires pthread")
AP.add_argument("--lrt", required=False, action='store_true', help="Requires lrt")
AP.add_argument("--nonative", required=False, action='store_true', help="Don't use -march=native (for OSX M1)")
AP.add_argument("--symbols", required=False, action='store_true', help="Debug symbols (default if --debug)")
AP.add_argument("--nostrip", required=False, action='store_true', help="Don't strip symbols (default if --debug or --profile or --symbols)")
AP.add_argument("--nostatic", required=False, action='store_true', help="Don't do static linking")
AP.add_argument("--nomake", required=False, action='store_true', help="Generate Makefile/make.bash only, don't run make")
AP.add_argument("--bash", required=False, action='store_true', help="Generate make.bash (default Makefile)")
AP.add_argument("--git_hash", required=False, action='store_true', help="Generate git_hash.h (default gitver.txt)")
AP.add_argument("--ec2", required=False, action='store_true', help="EC2 compatible -march=x86-64-v3 -mtune=native")
AP.add_argument("--x86-64-v2", required=False, action='store_true', help="-march=x86-64-v2, link d5 build-wsl-v2 (default)")
AP.add_argument("--x86-64-v3", required=False, action='store_true', help="-march=x86-64-v3, link d5 build-wsl-v3")
AP.add_argument("--native", required=False, action='store_true', help="-march=native, link d5 build-wsl-v3")

Args = AP.parse_args()
debug = Args.debug
profile = Args.profile
nomake = Args.nomake
std = Args.std
cppcompiler = Args.cppcompiler
ccompiler = Args.ccompiler
nostrip = profile or debug or Args.symbols
symbols = profile or debug or Args.symbols
static = True
if Args.nostatic or Args.santhread or Args.sanaddr:
	static = False
	symbols = True

deletes = []
if not Args.deletes is None:
	deletes = Args.deletes.split('+')

ProjFileName = None
HdrNames = []
for FileName in os.listdir("."):
	if FileName.endswith(".vcxproj"):
		ProjFileName = FileName
	elif FileName.endswith(".h"):
		HdrNames.append(FileName)
if ProjFileName is None:
	sys.stderr.write("\nProject file not found in current directory\n")
	sys.exit(1)

if Args.binary is None:
	binary = ProjFileName.replace(".vcxproj", "")
else:
	binary = Args.binary
	if binary.find('/') >= 0:
		sys.stderr.write("ERROR path name not allowed for --binary")
		sys.exit(0)
sys.stderr.write("binary=" + binary + "\n")

compiler_opts = " -iquote . -I%s -flto -ffast-math" % D5_MYUTILS_INC
linker_opts = " -flto -ffast-math"

if Args.x86_64_v2 and Args.x86_64_v3:
	sys.stderr.write("ERROR --x86-64-v2 and --x86-64-v3 are mutually exclusive\n")
	sys.exit(1)

# Portable Linux default is v2. Muscle has no AVX2 kernel; v3 gain is small.
if Args.ec2:
	D5_MYUTILS_LIB = D5_MYUTILS_LIB_V3
	compiler_opts += " -march=x86-64-v3 -mtune=native"
	linker_opts += " -march=x86-64-v3 -mtune=native"
elif Args.x86_64_v3:
	D5_MYUTILS_LIB = D5_MYUTILS_LIB_V3
	compiler_opts += " -march=x86-64-v3"
	linker_opts += " -march=x86-64-v3"
elif Args.native:
	D5_MYUTILS_LIB = D5_MYUTILS_LIB_V3
	compiler_opts += " -march=native"
	linker_opts += " -march=native"
elif Args.nonative:
	D5_MYUTILS_LIB = D5_MYUTILS_LIB_V2
else:
	D5_MYUTILS_LIB = D5_MYUTILS_LIB_V2
	compiler_opts += " -march=x86-64-v2"
	linker_opts += " -march=x86-64-v2"

sys.stderr.write("d5_lib=" + D5_MYUTILS_LIB + "\n")

if Args.santhread:
	compiler_opts += " -fsanitize=thread"
	linker_opts += " -fsanitize=thread"

if Args.sanaddr:
	compiler_opts += " -fsanitize=address"
	linker_opts += " -fsanitize=address"

if Args.profile:
	compiler_opts += " -fno-omit-frame-pointer"
	linker_opts += " -fno-omit-frame-pointer"

if std:
	compiler_opts += " --std=" + std

if debug:
	compiler_opts += " -O0 -DDEBUG"
	linker_opts += " -O0"
else:
	compiler_opts += " -O3 -DNDEBUG"
	linker_opts += " -O3"

if symbols:
	compiler_opts += " -g3"
	linker_opts += " -g3"

if Args.openmp:
	compiler_opts += " -fopenmp"
	linker_opts += " -fopenmp"

if Args.pthread:
	compiler_opts += " -pthread"
	linker_opts += " -lpthread"

rc = os.system('test -z $(git status --porcelain) 2> /dev/null')
if rc != 0:
	sys.stderr.write("\n\nWarning -- Uncommited changes\n\n")

rc = os.system(r'echo \"$(git log --oneline | head -n1 | cut "-d " -f1)\" | tee gitver.txt')
if rc != 0:
	sys.stderr.write("\n\nERROR -- failed to generate gitver.txt\n\n")
	sys.exit(1)
sys.stderr.write("gitver.txt done.\n")
if Args.git_hash:
	rc = os.system(r'h=`cat gitver.txt`; echo "#define GIT_HASH $h" > git_hash.h ; rm -f gitver.txt')
	if rc != 0:
		sys.stderr.write("\n\nERROR -- failed to generate gitver.txt\n\n")
		sys.exit(1)

rc = os.system(r'rm -rf o/ ../bin/%s*' % binary)
if rc != 0:
	sys.stderr.write("\n\nERROR -- failed to clean\n\n")
	sys.exit(1)
sys.stderr.write("clean done.\n")

OBJDIR = "o"
BINDIR = "../bin"


CPPNames = []
CNames = []
with open(ProjFileName) as File:
	for Line in File:
		Line = Line.strip()
		Line = Line.replace('"', '')
		Line = Line.replace(' ', '')
		if Line.startswith("<ClCompileInclude"):
			Fields = Line.split("=")
			if len(Fields) != 2:
				continue
			FileName = Fields[1]
			FileName = FileName.replace("/>", "")
			if FileName in deletes:
				sys.stderr.write("Deleting %s\n" % FileName)
				continue
			if FileName.endswith(".cpp"):
				FileName = FileName.replace(".cpp", "")
				CPPNames.append(FileName)
			elif FileName.endswith(".c"):
				FileName = FileName.replace(".c", "")
				CNames.append(FileName)

assert len(CPPNames) > 0 or len(CNames) > 0

if Args.bash:
	bashfn = "./make.bash.tmp"
	with open(bashfn, "w") as f:
		f.write("#/bin/bash -e\n")
		f.write("rm -rf o/\n")
		f.write("mkdir -p o/\n")
		n = len(CNames)
		if n > 0:
			f.write("for cname in \\\n")
			for i in range(n):
				f.write("   " + CNames[i])
				if i+1 < n:
					f.write(" \\");
				f.write("\n")
			f.write("do\n")
			f.write("   echo $cname >> o/cnames.tmp\n");
			f.write("done\n")
		n = len(CPPNames)
		if n > 0:
			f.write("for cppname in \\\n")
			for i in range(n):
				f.write("   " + CPPNames[i])
				if i+1 < n:
					f.write(" \\");
				f.write("\n")
			f.write("do\n")
			f.write("   echo $cppname >> o/cppnames.tmp\n");
			f.write("done\n")

		f.write("echo '#!/bin/bash -x' > o/compile_c.bash\n")
		f.write("echo '%s -c %s $1.c -o o/$1.o' >> o/compile_c.bash\n" % (ccompiler,  compiler_opts))
		f.write("echo '#!/bin/bash -x' > o/compile_cpp.bash\n")
		f.write("echo '%s -c %s $1.cpp -o o/$1.o' >> o/compile_cpp.bash\n" % (cppcompiler,  compiler_opts))
		f.write("chmod +x o/compile_c.bash  o/compile_cpp.bash\n")
		f.write("cat o/cnames.tmp \\\n")
		f.write("	| parallel o/compile_c.bash\n")
		f.write("cat o/cppnames.tmp \\\n")
		f.write("	| parallel o/compile_cpp.bash\n")
		f.write("%s %s \\\n" % (cppcompiler, linker_opts))
		for CName in CNames:
			f.write("	o/%s.o \\\n" % CName)
		for CPPName in CPPNames:
			f.write("	o/%s.o \\\n" % CPPName)
		f.write("	%s \\\n" % D5_MYUTILS_LIB)
		if static:
			f.write("	-static \\\n")
		f.write("	-o ../bin/%s\n" % binary)
		f.write("ls -lh ../bin/%s\n" % binary)
	os.system("chmod +x " + bashfn)
	rc = 0
	if not Args.nomake:
		sys.stderr.write("Running make...\n")
		rc = os.system(bashfn)
		sys.stderr.write("make done rc=%d\n" % rc)
	sys.stderr.write("exit rc=%d\n" % rc)
	sys.exit(rc)

with open("Makefile", "w") as f:
	def Out(s):
		f.write(s + "\n")

	BINPATH = "$(BINDIR)/%s" % (binary) 

	Out("######################################################")
	Out("# Makefile is generated by " + sys.argv[0])
	Out("# Don't edit the Makefile -- update the python script")
	Out("######################################################")
	Out("")
	Out("BINDIR := %s" % BINDIR)
	Out("OBJDIR := %s" % OBJDIR)
	Out("BINPATH := %s" % BINPATH)

	if CNames:
		Out("")
		Out("CC = " + ccompiler)
		Out("CFLAGS := " + compiler_opts)

	if CPPNames:
		Out("")
		Out("CPP = " + cppcompiler)
		Out("CPPFLAGS := " + compiler_opts)

	Out("")
	Out("UNAME_S := $(shell uname -s)")
	Out("LDFLAGS := $(LDFLAGS) " + linker_opts)
	if static:
		Out("ifeq ($(UNAME_S),Linux)")
		Out("    LDFLAGS += -static")
		Out("endif")

	Out("")
	Out("HDRS = \\")
	for Name in sorted(HdrNames):
		Out("  %s \\" % Name)

	Out("")
	Out("OBJS = \\")
	for Name in CPPNames:
		Out("  $(OBJDIR)/%s.o \\" % (Name))

	for Name in CNames:
		Out("  $(OBJDIR)/%s.o \\" % (Name))

	Out("")
	Out(".PHONY: clean")

	Out("")
	Out("$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)")

	if len(CPPNames) > 0:
		Cmd = "\t$(CPP) $(LDFLAGS) $(OBJS) %s -o $(BINPATH)" % D5_MYUTILS_LIB
	else:
		Cmd = "\t$(CC) $(LDFLAGS) $(OBJS) %s -o $(BINPATH)" % D5_MYUTILS_LIB

	if Args.lrt:
		Cmd += " -lrt"
	Out(Cmd)

	if not nostrip:
		Out("	strip $(BINPATH)")

	Out("")
	Out("$(OBJDIR)/ :")
	Out("	mkdir -p $(OBJDIR)/")

	Out("")
	Out("$(BINDIR)/ :")
	Out("	mkdir -p $(BINDIR)/")

	if CNames:
		Out("")
		Out("$(OBJDIR)/%.o : %.c $(HDRS)")
		Out("	$(CC) $(CFLAGS) -c -o $@ $<")

	if CPPNames:
		Out("")
		Out("$(OBJDIR)/%.o : %.cpp $(HDRS)")
		Out("	$(CPP) $(CPPFLAGS) -c -o $@ $<")

sys.stderr.write("Makefile done.\n")
if nomake:
	sys.exit(0)

sys.stderr.write("Running make...\n")
rc = os.system("make -j64 2> make.stderr | tee make.stdout")
sys.stderr.write("make done rc=%d\n" % rc)
os.system("tail make.stderr")
if rc != 0:
	sys.stderr.write("\n\nERROR -- make failed, see make.stderr\n\n")
	sys.exit(1)
sys.stderr.write("make done.\n")
os.system("ls -lh ../bin/" + binary + "\n")
