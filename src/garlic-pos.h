#ifndef __GARLIC_POS_H__
#define __GARLIC_POS_H__

#include <cstdint>

//Type of a physical (base pair) coordinate.
//
//This was int, so nothing past 2,147,483,647 bp could be represented.  That
//limit is reachable: several plant and amphibian chromosomes exceed it (Paris
//japonica, Pinus, Ambystoma), as does any workflow using concatenated or
//lifted coordinates.
//
//It did not fail loudly.  Positions are parsed with strtod and were then
//assigned to an int, and converting an out-of-range double to an integer type
//is undefined behaviour -- so the result is whatever the toolchain does.
//Measured here (clang, arm64) on chr21 shifted by 3,000,000,000 bp: every
//position saturated at INT_MAX, so all 171 ROH were emitted with the SAME
//coordinates (chromStart 2,147,483,646, chromEnd 2,147,483,647, length 1) and
//the program exited 0.  A different toolchain may wrap negative instead.
//With pos_t the same input gives 158 distinct starts spanning
//3,014,137,684..3,046,847,315 and chromEnd - chromStart == length throughout.
//
//It lives in its own header because garlic-data.h and garlic-centromeres.h
//both need it and garlic-data.h already includes garlic-centromeres.h.
typedef int64_t pos_t;

#endif
