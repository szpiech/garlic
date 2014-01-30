# Makefile for rohscan

#G++FLAG = -g -Wall
G++FLAG = -O3 -m64 -mmmx -msse -msse2
I_PATH = ../include
L_PATH = ../lib
#L_PATH2 = /usr/local/lib

#For static linking
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a $(L_PATH)/libann_figtree_version.a $(L_PATH)/libfigtree.a


all : rohscan

rohscan : rohscan-main.o rohscan-data.o param_t.o rohscan-roh.o rohscan-kde.o
	g++ -o rohscan rohscan-main.o rohscan-data.o param_t.o rohscan-roh.o rohscan-kde.o -pthread $(LINK_OPTS) -L$(L_PATH)

rohscan-main.o : rohscan-main.cpp
	g++ $(G++FLAG) -c rohscan-main.cpp -I$(I_PATH)

rohscan-data.o : rohscan-data.cpp
	g++ $(G++FLAG) -c rohscan-data.cpp 

rohscan-roh.o : rohscan-roh.cpp
	g++ $(G++FLAG) -c rohscan-roh.cpp 

rohscan-kde.o : rohscan-kde.cpp
	g++ $(G++FLAG) -c rohscan-kde.cpp -I$(I_PATH)

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

clean :
	rm *.o
