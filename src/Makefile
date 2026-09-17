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
  # MSYS2/MinGW report MINGW64_NT-10.0-..., MSYS_NT-... or CYGWIN_NT-...  None
  # of those is Darwin, so they used to fall through to PLATFORM=linux and the
  # build then looked for ../lib/linux and passed -static -Wl,--whole-archive.
  # The manual says "type make; the Makefile detects the platform from uname,
  # so no editing is needed" -- that was only true off Windows.
  else ifneq (,$(findstring MINGW,$(UNAME_S)))
    PLATFORM = win32
  else ifneq (,$(findstring MSYS,$(UNAME_S)))
    PLATFORM = win32
  else ifneq (,$(findstring CYGWIN,$(UNAME_S)))
    PLATFORM = win32
  else
    PLATFORM = linux
  endif
endif

# -- optimisation ------------------------------------------------------------
# Portable and safe by default.  Add -march=native yourself for a local build.
OPTFLAGS ?= -O3
WARNFLAGS ?= -Wall
STDFLAG ?= -std=c++11

# Stamped into --version and into the git_sha field of <out>.params.json.
# "unknown" when built from a source tarball with no .git.
#
# This used to be passed as -D on the garlic-cli.o rule, which meant the stamp
# only refreshed when garlic-cli.cpp itself changed: --version and the
# provenance record in params.json both reported whatever commit was current
# the last time that one file was touched.  It now goes through a generated
# header that garlic-cli.o depends on, rewritten only when the value actually
# changes so it does not force a rebuild on every make.
GIT_SHA := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)

# -- per-platform settings ---------------------------------------------------
# There is no L_PATH any more.  Every platform now takes zlib from the
# toolchain, so nothing is linked out of lib/ and the -L that used to be on the
# link line is gone.  It was not merely redundant: it put ../lib/win32 on the
# search path, where the vendored libz.a is a 32-bit i386 archive (COFF machine
# 0x014c, as is libpthreadGC2.a and the garlic.exe under bin/win32).  A 64-bit
# MinGW then skipped it as incompatible and the link failed on undefined zlib
# symbols -- which is exactly how the first Windows CI run failed.

# GSL is no longer used.  garlic-math.cpp provides the ten routines it was
# linked for, reproducing each bit for bit (see that file's header), so the
# per-platform libgsl.a/libgslcblas.a archives are not needed and neither is
# a system GSL.  The only remaining external dependency is zlib -- from the
# system on Linux and macOS, from lib/win32/libz.a on win32.  Threading is
# std::thread, so no pthreads library is linked by name on any platform
# (-pthread on linux is a toolchain flag, not a vendored archive).
# figtree/ANN are no longer linked: the KDE evaluates its Gaussian sum directly
# (see kdeGaussian in garlic-kde.cpp).  The archives remain under lib/ for the
# moment; nothing references them.

ifeq ($(PLATFORM),linux)
  ARCHFLAGS = -m64 -pthread
  # -static is right for the bundled archives (the point of vendoring them is a
  # binary that runs on a cluster without matching libraries).  Build with
  # STATIC= to turn it off.
  STATIC ?= -static
  #
  # std::thread needs more care under -static than pthread_create did.
  # libstdc++ decides at RUNTIME whether threading is available by testing a
  # WEAK reference (__gthread_active_p, historically via pthread_key_create).
  # In a static link the linker only pulls in the archive members something
  # references strongly, so that symbol can come out unresolved-but-weak, the
  # test reports single-threaded, and the first std::thread constructor throws
  #   std::system_error: Enable multithreading to use std::thread
  # at runtime even though the link succeeded.  Direct pthread_create calls
  # did not have this problem, because calling it IS a strong reference.
  # --whole-archive forces all of libpthread in and restores the check.
  #
  # On glibc 2.34+ libpthread is merged into libc and this is moot, but -static
  # exists here for old cluster systems, which is exactly where it is not.
  ifeq ($(STATIC),-static)
    LINK_OPTS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
  else
    LINK_OPTS = $(STATIC) -pthread
  endif
else ifeq ($(PLATFORM),osx)
  # x86_64 macOS: SSE2 is baseline on x86_64, kept explicit for older toolchains
  ARCHFLAGS = -m64 -msse2
  LINK_OPTS =
else ifeq ($(PLATFORM),macos-arm)
  # arm64: the x86 ISA flags (-m64 -mmmx -msse -msse2) are meaningless here
  ARCHFLAGS =
  LINK_OPTS =
else ifeq ($(PLATFORM),win32)
  # No -I ../include/win32 either: that directory holds the 2017 zlib headers
  # that went with the 32-bit archive.  MSYS2 supplies zlib.h and libz for the
  # architecture actually being built (mingw-w64-x86_64-zlib for MINGW64,
  # mingw-w64-i686-zlib for MINGW32), which is the only combination that can
  # work -- a 64-bit build cannot use an i386 archive at all.
  # -DPTW32_STATIC_LIB is gone: it guarded the pthread_win32_process_attach_np
  # /detach_np pair in main, which pthreads-win32 requires and std::thread does
  # not.  Nothing tests the macro now.
  #
  # lib/win32/libpthreadGC2.a is also gone from the link.  That archive is
  # pthreads-win32, a separate project from the winpthreads that MinGW-w64's
  # libstdc++ builds std::thread on, so it is not what a std::thread binary
  # needs: a MinGW build now requires a GCC configured with the POSIX threading
  # model, which supplies threading itself.  With a win32-threads GCC,
  # libstdc++ does not define std::thread at all and the build fails at compile
  # time rather than producing a subtly single-threaded binary.
  #
  # UNTESTED: no Windows toolchain is available to me, and CI has no Windows
  # runner.  zlib is still linked from lib/win32.
  ARCHFLAGS =
  # zlib comes from -lz on the link line like every other platform now.
  LINK_OPTS = -static-libgcc -static-libstdc++
else
  $(error Unknown PLATFORM '$(PLATFORM)'; expected linux, osx, macos-arm or win32)
endif

CXXFLAGS ?= $(OPTFLAGS) $(WARNFLAGS) $(STDFLAG) $(ARCHFLAGS)
G++FLAG = $(CXXFLAGS)
CC = $(CXX)

# MinGW's ld appends .exe to an output name that has none, so `-o garlic`
# produces garlic.exe and the target `garlic` is never created.  Make then
# thinks the target is missing and relinks on every invocation, and `clean`
# leaves the binary behind.  It only appeared to work because the MSYS2
# runtime resolves a stat of "garlic" to "garlic.exe"; naming the target
# correctly means make's dependency tracking is right rather than accidentally
# right.
ifeq ($(PLATFORM),win32)
  EXE = .exe
else
  EXE =
endif


all : garlic$(EXE)

OBJS = garlic-main.o garlic-data.o param_t.o garlic-roh.o garlic-kde.o \
       garlic-cli.o garlic-errlog.o gmm.o BoundFinder.o gzstream.o \
       garlic-centromeres.o garlic-pbar.o garlic-math.o \
       garlic-options.o

garlic$(EXE) : $(OBJS)
	$(CC) -o garlic$(EXE) $(OBJS) -lz $(LINK_OPTS)

# Header dependencies are generated by the compiler (-MMD -MP) rather than
# written out by hand.  The hand-written lists were incomplete -- garlic-main.o
# was declared to depend only on garlic-cli.h, though garlic-main.cpp includes
# garlic-data.h, garlic-roh.h, garlic-kde.h, param_t.h and
# garlic-centromeres.h.  Editing a struct in garlic-data.h therefore rebuilt
# NOTHING, and the link silently produced a binary whose translation units
# disagreed about that struct's layout: adding a field to IndData gave a
# segfault that went away on 'make clean'.  -MMD writes a .d file per object
# listing every header it actually read; -MP adds phony targets so a deleted
# header does not break the build.
DEPFLAGS = -MMD -MP

%.o : %.cpp
	$(CC) $(G++FLAG) $(DEPFLAGS) -c $< -I$(I_PATH) $(I_PATH2)

# garlic-version.h is generated, so it has to exist before the first compile
# that includes it; the .d file cannot know about it on a clean tree.
garlic-cli.o : garlic-version.h

# Written on every make; replaced only when the SHA changes.
garlic-version.h : FORCE
	@printf '/* generated by the Makefile -- do not edit */\n#define GARLIC_GIT_SHA "%s"\n' '$(GIT_SHA)' > $@.tmp
	@cmp -s $@.tmp $@ 2>/dev/null && rm -f $@.tmp || mv $@.tmp $@

FORCE :

-include $(OBJS:.o=.d)

clean :
	rm -f *.o *.d garlic garlic.exe garlic-version.h garlic-version.h.tmp


# Regenerate the command line reference in README and manual/flags.tex from the
# program's own help strings, so the documentation cannot drift from the binary.
# The flag the manual described that garlic never accepted (--auto-winsize-size,
# whose body actually documented --auto-winsize-step), and the nineteen the
# program gained without being documented at all, were both consequences of
# maintaining that list by hand.
docs : garlic
	./garlic --dump-docs tex > ../manual/flags.tex
	awk '/^#BEGIN GENERATED FLAG REFERENCE/{print; system("./garlic --dump-docs txt"); s=1; next} \
	     /^#END GENERATED FLAG REFERENCE/{s=0} !s' ../README > README.docs.tmp
	mv README.docs.tmp ../README
	@echo "regenerated ../README flag block and ../manual/flags.tex"

.PHONY : all clean docs test FORCE

# Run the test suite: unit tests, golden-output regression on the bundled
# example, determinism, the exit-code contract, the --load-params round trip,
# and BED conformance.  Needs only a shell and a C++ compiler.
test : garlic
	sh ../test/run_tests.sh
