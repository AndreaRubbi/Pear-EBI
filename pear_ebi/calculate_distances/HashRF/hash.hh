//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
//
//		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		HashRFMap: Class for hashing bitstrings
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

#ifndef HASHRFMYHASH_HH
#define HASHRFMYHASH_HH

#include <map>
#include <vector>

#include "hashfunc.hh"

using namespace std;

typedef struct {
	unsigned long long  _hv2;
	unsigned int 	      _t_i;
	float 				      _dist;
} BUCKET_STRUCT_T;

typedef vector<BUCKET_STRUCT_T> V_BUCKET_T;
typedef vector<V_BUCKET_T> HASHTAB_T;


typedef struct {
	unsigned long	long		_hv2;
	vector<unsigned int> 	_vec_treeidx;
	vector<float>     		_vec_dist;
} TREEIDX_STRUCT_T;

typedef vector<TREEIDX_STRUCT_T> V_BUCKET_T2;
typedef vector<V_BUCKET_T2> HASHTAB_T2;


class HashRFMap {

public:
	HashRFMap() {}
	~HashRFMap() {}

	CHashFunc 	_HF;
	HASHTAB_T 	_hashtab;	// for weighted RF
	HASHTAB_T2 	_hashtab2;

	void hashing_bs_without_type2_nbits(unsigned int tree_i, unsigned int num_taxa, unsigned long long hv1, unsigned long long hv2, float dist, bool w_option); // fast hash-rf
	void uhashfunc_init(unsigned int t, unsigned int n, unsigned int c);

	void hashrfmap_clear();
};


#endif
