# Makefile for garlic

CC = g++
#CC = clang++
#G++FLAG = -g -Wall

I_PATH = ../include

#osx
L_PATH = ../lib/osx
G++FLAG = -fast -m64 -mmmx -msse -msse2
#For static linking
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/libann_figtree_version.a $(L_PATH)/libfigtree.a

#linux
#L_PATH = ../lib/linux
#G++FLAG = -O3 -m64 -mmmx -msse -msse2
#For static linking on linux boxes
#LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/libann_figtree_version.a $(L_PATH)/libfigtree.a -static -pthread

#win
#I_PATH2 = -I ../include/win32
#L_PATH = ../lib/win32
#G++FLAG = -O3 -DPTW32_STATIC_LIB
#For static linking
#LINK_OPTS = $(L_PATH)/libpthreadGC2.a $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/ann_figtree_version.lib $(L_PATH)/figtree.lib $(L_PATH)/libz.a -static-libgcc -static-libstdc++


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
	rm *.o
