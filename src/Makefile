# The $(OS) variable is the o/s name returned by uname, which is 
# used as the sub-directory name under src/ where object files
# and the executable are stored. This allows several target
# operating systems in the same directory structure.
# Typical values are:

#	Platform    Value of $(OS)
#   --------    --------------
#	Linux       linux
#	Mac OSX     Darwin
#	Cygwin      CYGWIN_NT-10.0

# Building on Mac OSX is challenging because Apple does not support gcc or
# the OMP threading library. Hacks to install gcc and OMP vary by OSX release.
# This Makefile works with the AWS Catalina v10.15.7 AMI. With this AMI, 
# running 'brew install gcc' currently installs gcc v11.

OS := $(shell uname)

CPPFLAGS += -DNDEBUG -pthread

# Detect if CXX is Apple clang, use custom preprocessor call
COMPILER_VERSION := $(shell $(CXX) --version)
ifneq '' '$(findstring Apple clang, $(COMPILER_VERSION))'
	CPPFLAGS += -Xpreprocessor -fopenmp
else
	CPPFLAGS += -fopenmp
endif

# If Darwin, then -fopenmp won't work in linker stage appropriately.
# (Apple-ld is picky.)
# As such, the user must specify an OpenMP implementation in LIBS.
# We check this (quick-and-dirty) to help them out.
ifneq '' '$(findstring Darwin, $(OS))'
	ifeq '' '$(findstring omp, $(LIBS))'
		maybe_missing_omp=1
	endif
else
	LDFLAGS += -fopenmp
endif

# Add default flags
CXXFLAGS += -O3 -ffast-math -std=c++11
LDFLAGS += -O3 -pthread -lpthread

HDRS := $(shell echo *.h)
OBJS := $(shell echo *.cpp | sed "-es/^/$(OS)\//" | sed "-es/ / $(OS)\//g" | sed "-es/\.cpp/.o/g")
SRCS := $(shell ls *.cpp *.h)


.PHONY: clean

$(OS)/muscle : gitver.txt $(OS)/ $(OBJS)
	@if [ -n "${maybe_missing_omp}" ]; then \
		echo "Warning: please specify an OpenMP library in LIBS."; \
	fi

	$(CXX) $(LDFLAGS) $(LIBS) $(OBJS) -o $@

	@# Warning: do not add -d option to strip, this is not portable
	strip $(OS)/muscle

gitver.txt : $(SRCS)
	bash ./gitver.bash

$(OS)/ :
	mkdir -p $(OS)/

$(OS)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf gitver.txt $(OS)/
