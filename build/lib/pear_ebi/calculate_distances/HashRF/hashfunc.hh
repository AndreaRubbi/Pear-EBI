//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		CHashFunc: Universal hash functions
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).
//
//*****************************************************************/

#ifndef HASHFUNC_HH
#define HASHFUNC_HH

#include <iostream>
#include <fstream>


typedef struct {
	unsigned long long hv1;
	unsigned long long hv2;
} HV_STRUCT_T;

class CHashFunc {

	unsigned long long    _m1; // prime number1 for hash function1
	unsigned long long    _m2; // prime number1 for hash function2
	unsigned int          _t;  // number of trees
	unsigned int          _n;  // number of taxa
	unsigned long long *  _a1; // random numbers for hash function1
	unsigned long long *  _a2; // random numbers for hash function2
	unsigned int		      _c;	 // double collision factor: constant for c*t*n of hash function2;

public:
	CHashFunc() : _m1(0), _m2(0), _t(0), _n(0), _a1(NULL), _a2(NULL), _c(0) {}
	CHashFunc(unsigned int t, unsigned int n, unsigned int c);
  ~CHashFunc();

	void UHashfunc_init(unsigned int t, unsigned int n, unsigned int c);
	void UHashFunc(HV_STRUCT_T &hv, uint64_t bs64, unsigned numBits);

	unsigned long long GetPrime(unsigned long long topNum, unsigned from);

	// Implicit bp
	unsigned long long getA1(unsigned idx) { return (_a1[idx]); }
	unsigned long long getA2(unsigned idx) { return (_a2[idx]); }
	unsigned long long getM1() { return _m1; }
	unsigned long long getM2() { return _m2; }
};

#endif
