
#AM_CXXFLAGS = -static
AM_CXXFLAGS = -Wall
#AM_CXXFLAGS += -pg
AM_CXXFLAGS += -O2 -funroll-loops -fexpensive-optimizations -DNDEBUG

# real program stuff
bin_PROGRAMS = hashrf

hashrf_SOURCES =
hashrf_SOURCES += hashfunc.hh hashfunc.cc
hashrf_SOURCES += label-map.hh label-map.cc
hashrf_SOURCES += hash.hh hash.cc

### Newick parser
hashrf_SOURCES += newick.h newick.c

# Mersenne Twister
hashrf_SOURCES += Random.cpp Random.hpp


hashrf_SOURCES += hashrf.cc

#eof
