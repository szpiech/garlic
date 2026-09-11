# Makefile for garlic
#
# Platform (library path and link options) is detected automatically from
# `uname`.  Override any variable on the command line, e.g.
#
#     make                       # autodetect
#     make CXX=clang++           # different compiler
#     make OPTFLAGS='-O2 -g'     # different optimisation level
#     make PLATFORM=win32        # cross/explicit platform
#
# NOTE: do not use `-fast` here.  It is an Intel/legacy-Apple spelling that
# clang accepts and silently ignores ("argument unused during compilation"),
# which left the macOS build at -O0 and ~3.4x slower than -O2.

CXX ?= g++
I_PATH = ../include

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifndef PLATFORM
  ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
      PLATFORM = macos-arm
    else
      PLATFORM = osx
    endif
  else
    PLATFORM = linux
  endif
endif

# -- optimisation ------------------------------------------------------------
# Portable and safe by default.  Add -march=native yourself for a local build.
OPTFLAGS ?= -O3
WARNFLAGS ?= -Wall
STDFLAG ?= -std=c++11

# -- per-platform settings ---------------------------------------------------
L_PATH = ../lib/$(PLATFORM)
GSL_LIBS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a
# figtree/ANN are no longer linked: the KDE evaluates its Gaussian sum directly
# (see kdeGaussian in garlic-kde.cpp).  The archives remain under lib/ for the
# moment; nothing references them.

ifeq ($(PLATFORM),linux)
  ARCHFLAGS = -m64
  LINK_OPTS = $(GSL_LIBS) -static -pthread
else ifeq ($(PLATFORM),osx)
  # x86_64 macOS: SSE2 is baseline on x86_64, kept explicit for older toolchains
  ARCHFLAGS = -m64 -msse2
  LINK_OPTS = $(GSL_LIBS)
else ifeq ($(PLATFORM),macos-arm)
  # arm64: the x86 ISA flags (-m64 -mmmx -msse -msse2) are meaningless here
  ARCHFLAGS =
  LINK_OPTS = $(GSL_LIBS)
else ifeq ($(PLATFORM),win32)
  I_PATH2 = -I ../include/win32
  ARCHFLAGS = -DPTW32_STATIC_LIB
  LINK_OPTS = $(L_PATH)/libpthreadGC2.a $(GSL_LIBS) \
              $(L_PATH)/libz.a -static-libgcc -static-libstdc++
else
  $(error Unknown PLATFORM '$(PLATFORM)'; expected linux, osx, macos-arm or win32)
endif

CXXFLAGS ?= $(OPTFLAGS) $(WARNFLAGS) $(STDFLAG) $(ARCHFLAGS)
G++FLAG = $(CXXFLAGS)
CC = $(CXX)


all : garlic

garlic : garlic-main.o garlic-data.o param_t.o garlic-roh.o garlic-kde.o garlic-cli.o garlic-errlog.o gmm.o BoundFinder.o gzstream.o garlic-centromeres.o garlic-pbar.o
	$(CC) -o garlic garlic-main.o garlic-data.o param_t.o garlic-roh.o garlic-kde.o garlic-cli.o garlic-errlog.o gmm.o BoundFinder.o gzstream.o garlic-centromeres.o garlic-pbar.o -lz $(LINK_OPTS) -L$(L_PATH)

garlic-main.o : garlic-main.cpp garlic-cli.h
	$(CC) $(G++FLAG) -c garlic-main.cpp -I$(I_PATH) $(I_PATH2)

garlic-data.o : garlic-data.cpp garlic-data.h
	$(CC) $(G++FLAG) -c garlic-data.cpp  -I$(I_PATH) $(I_PATH2)

garlic-roh.o : garlic-roh.cpp garlic-roh.h
	$(CC) $(G++FLAG) -c garlic-roh.cpp -I$(I_PATH) $(I_PATH2)

garlic-kde.o : garlic-kde.cpp garlic-kde.h
	$(CC) $(G++FLAG) -c garlic-kde.cpp -I$(I_PATH) $(I_PATH2)

garlic-cli.o : garlic-cli.cpp garlic-cli.h
	$(CC) $(G++FLAG) -c garlic-cli.cpp -I$(I_PATH) $(I_PATH2)

garlic-errlog.o : garlic-errlog.cpp garlic-errlog.h
	$(CC) $(G++FLAG) -c garlic-errlog.cpp -I$(I_PATH) $(I_PATH2)

garlic-centromeres.o : garlic-centromeres.cpp garlic-centromeres.h
	$(CC) $(G++FLAG) -c garlic-centromeres.cpp -I$(I_PATH) $(I_PATH2)

garlic-pbar.o : garlic-pbar.cpp garlic-pbar.h
	$(CC) $(G++FLAG) -c garlic-pbar.cpp -I$(I_PATH) $(I_PATH2)

param_t.o : param_t.cpp param_t.h
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp gzstream.h
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH2)

gmm.o : gmm.cpp gmm.h
	$(CC) $(G++FLAG) -c gmm.cpp -I$(I_PATH) $(I_PATH2)

BoundFinder.o : BoundFinder.cpp BoundFinder.h
	$(CC) $(G++FLAG) -c BoundFinder.cpp -I$(I_PATH) $(I_PATH2)

clean :
	rm -f *.o garlic

.PHONY : all clean
