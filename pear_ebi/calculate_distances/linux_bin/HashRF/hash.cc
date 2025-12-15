//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		HashRFMap: Class for hashing bit
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

#include "hash.hh"
#include "hashfunc.hh"


void
HashRFMap::uhashfunc_init(
	unsigned int t,
	unsigned int n,
//	float r,
	unsigned int c)
{
//	_HF.UHashfunc_init(t, n, r, c);
	_HF.UHashfunc_init(t, n, c);
}


void
HashRFMap::hashing_bs_without_type2_nbits(
	unsigned int treeIdx,
	unsigned int numTaxa,
	unsigned long long hv1,
	unsigned long long hv2,
	float dist,
	bool w_option)
{
	///////////////////////////////
	// double linked list
	///////////////////////////////
	unsigned sizeVec = _hashtab2[hv1].size();
	if (sizeVec > 0)
	{
		bool found = false;
		for (unsigned int i=0; i<sizeVec; ++i)
		{
			if (_hashtab2[hv1][i]._hv2 == hv2) {
				_hashtab2[hv1][i]._vec_treeidx.push_back(treeIdx);
				if (w_option)
					_hashtab2[hv1][i]._vec_dist.push_back(dist);

				found = true;
				break;
			}
		}
		if (!found)
		{
			TREEIDX_STRUCT_T bk2;
			bk2._hv2 = hv2;
			bk2._vec_treeidx.push_back(treeIdx);
			if (w_option)
				bk2._vec_dist.push_back(dist);
			_hashtab2[hv1].push_back(bk2);
		}
	}
	else if (sizeVec == 0)
	{
		TREEIDX_STRUCT_T bk2;
		bk2._hv2 = hv2;
		bk2._vec_treeidx.push_back(treeIdx);
		if (w_option)
			bk2._vec_dist.push_back(dist);
		_hashtab2[hv1].push_back(bk2);
	}
}

void
HashRFMap::hashrfmap_clear()
{
  _hashtab.clear();
	_hashtab2.clear();
}

// eof
