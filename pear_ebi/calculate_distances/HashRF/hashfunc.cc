//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		CHashFunc: Class for Universal hash functions
//
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

#include "hashfunc.hh"
#include <cassert>
#include "RandomLib/Random.hpp"

// constructor
CHashFunc::CHashFunc(
	unsigned int t,
	unsigned int n,
	unsigned int c)
	: _m1(0), _m2(0), _t(t), _n(n), _a1(NULL), _a2(NULL), _c(c)
{
	UHashfunc_init(t, n, c);
}

// destructor
CHashFunc::~CHashFunc()
{
	delete[] _a1;
	delete[] _a2;
}

void
CHashFunc::UHashfunc_init(
	unsigned int t,
	unsigned int n,
	unsigned int c)
{
	// Init member variables
	_t = t;
	_n = n;
	_c = c;

	// Get the prime number which is larger than t*n
	unsigned long long top = _t*_n;

	unsigned long long p = 0;
	unsigned int mul = 1;
	do {
		unsigned from = 100 * mul;
		p = GetPrime(top, from);
		++mul;
	} while (p == 0);

	_m1 = p;   	// t*n ~~ m1

	unsigned long long top2 = _c*_t*_n;
	unsigned long long p2 = 0;
	mul = 1;
	do {
		unsigned from = 100 * mul;
		p2 = GetPrime(top2, from);
		++mul;
	} while (p2 == 0);

	_m2 = p2; // m2 > c*t*n --> to avoid double collision ==> I just use _c*top for _m2
	_a1 = new unsigned long long[_n];
	_a2 = new unsigned long long[_n];


	// generate n random numbers between 0 and m1-1
	// for hash function1 and hash function2
	// rand() % 48
	// random number between 0~47
	RandomLib::Random rnd;		// r created with random seed

	for (unsigned int i=0; i<_n; ++i) {
		_a1[i] = rnd.Integer<unsigned long long>(_m1-1);
		_a2[i] = rnd.Integer<unsigned long long>(_m2-1);
	}
}



// Generate a prime number right after topNum
unsigned long long
CHashFunc::GetPrime(unsigned long long topNum, unsigned from)
{
	unsigned long long primeNum=0;
	unsigned long long candidate=0;

	if (topNum <= 100)
		candidate = 2;
	else
		candidate = topNum;

	while (candidate <= topNum+from) {
		unsigned long long trialDivisor = 2;
		int prime = 1;

		while (trialDivisor * trialDivisor <= candidate) {
			if (candidate % trialDivisor == 0) {
				prime = 0;
				break;
			}
			trialDivisor++;
		}
		if (prime) primeNum = candidate;

		candidate++;
	}

	return primeNum;
}
