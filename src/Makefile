# Makefile for rohscan

#G++FLAG = -g -Wall
G++FLAG = -O3 -m64 -mmmx -msse -msse2
I_PATH = ../include

#osx
L_PATH = ../lib/osx

#linux
#L_PATH = ../lib/linux

#For static linking
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/libann_figtree_version.a $(L_PATH)/libfigtree.a

#For static linking on linux boxes
#LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/libann_figtree_version.a $(L_PATH)/libfigtree.a -static


all : rohscan

rohscan : rohscan-main.o rohscan-data.o param_t.o rohscan-roh.o rohscan-kde.o gmm.o BoundFinder.o
	g++ -o rohscan rohscan-main.o rohscan-data.o param_t.o rohscan-roh.o rohscan-kde.o gmm.o BoundFinder.o -pthread $(LINK_OPTS) -L$(L_PATH)

rohscan-main.o : rohscan-main.cpp
	g++ $(G++FLAG) -c rohscan-main.cpp -I$(I_PATH)

rohscan-data.o : rohscan-data.cpp
	g++ $(G++FLAG) -c rohscan-data.cpp  -I$(I_PATH)

rohscan-roh.o : rohscan-roh.cpp
	g++ $(G++FLAG) -c rohscan-roh.cpp -I$(I_PATH)

rohscan-kde.o : rohscan-kde.cpp
	g++ $(G++FLAG) -c rohscan-kde.cpp -I$(I_PATH)

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

gmm.o : gmm.cpp
	g++ $(G++FLAG) -c gmm.cpp -I$(I_PATH)

BoundFinder.o : BoundFinder.cpp
	g++ $(G++FLAG) -c BoundFinder.cpp -I$(I_PATH)

clean :
	rm *.o
