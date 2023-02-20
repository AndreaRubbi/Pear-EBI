//*****************************************************************/
/*
		HashRF v.6.0.0

		<<<<<<<<<<<<<< 2007.10.19 >>>>>>>>>>>>>>>

		Fast algorithm to compute Robinson_Foulds topological distance between
		phylogenetic trees based on universal hashing.

	AUTHOR:
		Copyright (C) 2006, 2007 Seung-Jin Sul
			Dept. of Computer Science
			Texas A&M University
			U.S.A.
		(contact: sulsj@cs.tamu.edu)

  LICENSE AGREEMENT
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).

	DESCRIPTION:
		1. Read the first Newick tree and parsing it to tree data structure.
		2. Collect taxon names.
		3. Read all the trees and parsing them.
		4. Collect bipartitions and hash them.
		5. Compute RF distance.

	HISTORY:
		11/18/2006  Start to code new version
		11/19/2006  To read tree file and to parse trees, libcov library is used
		11/21/2006  Draw problem
		11/22/2006  First tree read problem == libcov covIO::ReadTree()

		11/30/2006  *CRITICAL MEM MGNT: dfs_traverse ==> delete bs


		12/01/2006  Add the constants, HASHTABLE_FACTOR and C
		12/02/2006  NUM_TREES and NUM_TAXA ==> make global for speedup
		12/03/2006  Reimplement GetTaxaLabels() and dfs_traverse()
		12/10/2006  Allocated memory clear

		01/01/2007  v 2.0.0 Remove opt library.
		01/02/2007  Remove STL Hash, jsut use my hash class

		01/17/2007  v.3.0.1 11000 + 00111 = 11111 --> no need to hash
					numBitstr --> count the number of bipartitions for each tree and
					check the number is less than (n-3).

		01/18/2007  delete tree one by one

		01/20/2007  std random number generator ==>  Mersenne Twister random
					number generator
					(http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html)

		02/01/2007  Print option #4 for ISMB07

		02/08/2007  After hashing a bs, the stack_bs should be cleared.
		02/08/2007  Reroot -> correct RF distance value

		03/17/2007  Incorporate the methods from libcov
		03/19/2007  Fix print option 1 & 2

		03/19/2007  Add the 64-bit representation scheme of PGM to Fast Hash-RF.
		03/19/2007  dfs_traverse_pgm() is implemented.

		03/22/2007  Some data files could not be read in Hash-RF.
		03/23/2007  Add myReadTree to covIO.cpp covIO.h
					To use covTree* instead of vector<covTree*>
					To debug corrupted double-linked list error

		03/26/2007  Four types of operations
					1. Without TYPE-III checking and n-bits representation  (Fast Hash-RF)
					2. Without TYPE-III checking and 64-bits representation (PGM+Fast Hash-RF)
					3. With    TYPE-III checking and n-bits representation  (Hash-RF)
					4. With    TYPE-III checking and 64-bits representation (PGM+Hash-RF)


		09.11.2007 Hash function init
					64bit --> map_hashrf.uhashfunc_init(NUM_TREES, BITS, HASHTABLE_FACTOR, C);
																												 ----

		09.11.2007 Hash function
					hv.hv1 += _a1[i]*ii;
					hv.hv2 += _a2[i]*ii;
					==>
					hv.hv1 += _a1[i];
					hv.hv2 += _a2[i];



		09.07.2007  64-bit optimization ????????????????????????????????????
					uint64_t ==> stack<uint64_t>
					vec_bucket = iter_hashcs->second; and vec_bucket[i]
					===> (iter_hashcs->second).size(); // remove copy



		09.25.2007 Multifurcating -> dissimilarity matrix ==> RF matrix

		10.18.2007 TCLAP argument processing



		10.22.2007 RF matrix init
					vector< vector<int> > DISSIM(NUM_TREES, vector<int>(NUM_TREES, NUM_TAXA-3));

		10.27.2007 Wrong RF values in 64bit modes.
					==> 64bit --> map_hashrf.uhashfunc_init(NUM_TREES, BITS, HASHTABLE_FACTOR, C);
					==> DO NOT USE

		10.30.2007 OPTIMIZATION
					Simple vector
					Implicit BP

		10.31.2007
					covTree <= DeleteAllNodes()
					Remove stacking



		11.12.2007 DEBUG
					In hashfunc.cc --> top2 = c*t*n;


		11.13.2007 data type
					int --> unsigned
					unsigned long

		11.14.2007 DEBUG
					DISSIM += ===> SIM -=

		11.28.2007 OPTIMIZATION
					Without bs.count_ones_zeros();
					No checking for nontrivial bipartition for implicit bipartition

		11.29.2007 OPTIMIZATION
					double linked list to optimize distance matrix computation.


		12.11.2007 Newick parser

		12.17.2007 Unweighted + weighted


		12.16.2008 Change matrix data type to unsigned short


		1.20.2009 ULLONG_MAX check

		2.13.2009 PROBLEM with c value and error rate
		    due to abnormal behavior of MT random number generator
		    fixed with RandomLib which is based on new MT random number generator

		2.17.2009 HashRF 6.0.0

		2.18.2009 remove bitset.hh bitset.cc

    4.14.2009 prepare distribution
        SIM -> divide into SHORTMAT & FLAOTMAT






*/
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "./tclap/CmdLine.h"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

#define LEFT               								0
#define RIGHT              								1
#define ROOT               								2
#define LABEL_WIDTH        								3

// Set a random number for m1 (= Initial size of hash table)
// m1 is the closest value to (t*n)+(t*n*HASHTABLE_FACTOR)
#define HASHTABLE_FACTOR    							0.2

// the c value of m2 > c*t*n in the paper
static unsigned int C						    			= 1000;
static unsigned NUM_TREES                 = 0; // number of trees
static unsigned NUM_TAXA                	= 0; // number of taxa
static bool WEIGHTED                      = false; // unweighted
static unsigned PRINT_OPTIONS							= 3; // matrix


#define __HALF_MAX_SIGNED(type) ((type)1 << (sizeof(type)*8-2))
#define __MAX_SIGNED(type) (__HALF_MAX_SIGNED(type) - 1 + __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type) (-1 - __MAX_SIGNED(type))

#define __MIN(type) ((type)-1 < 1?__MIN_SIGNED(type):(type)0)
#define __MAX(type) ((type)~__MIN(type))

#define assign(dest,src) ({ \
  typeof(src) __x=(src); \
  typeof(dest) __y=__x; \
  (__x==__y && ((__x<1) == (__y<1)) ? (void)((dest)=__y),0:1); \
})

#define add_of(c,a,b) ({ \
  typeof(a) __a=a; \
  typeof(b) __b=b; \
  (__b)<1 ? \
    ((__MIN(typeof(c))-(__b)<=(__a)) ? assign(c,__a+__b):1) : \
    ((__MAX(typeof(c))-(__b)>=(__a)) ? assign(c,__a+__b):1); \
})

void
GetTaxaLabels(
  NEWICKNODE *node,
  LabelMap &lm)
{
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  }
  else
    for (int i=0;i<node->Nchildren;++i)
      GetTaxaLabels(node->child[i], lm);
}

void
dfs_compute_hash(
  NEWICKNODE* startNode,
  LabelMap &lm,
  HashRFMap &vvec_hashrf,
  unsigned treeIdx,
  unsigned &numBitstr,
  unsigned long long m1,
  unsigned long long m2)
{
  // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
  // and push the bit string into stack
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = lm[temp];

    // Implicit BPs /////////////////////////
    // Set the hash values for each leaf node.
    startNode->hv1 = vvec_hashrf._HF.getA1(idx);
    startNode->hv2 = vvec_hashrf._HF.getA2(idx);
  }
  else {
    for (int i=0; i<startNode->Nchildren; ++i) {
      dfs_compute_hash(startNode->child[i], lm, vvec_hashrf, treeIdx, numBitstr, m1, m2);
    }

    // For weighted RF
    float dist = 0.0;
    if (WEIGHTED)
      dist = startNode->weight;
    else
      dist = 1;

    ++numBitstr;

    // Implicit BPs ////////////
    // After an internal node is found, compute the hv1 and hv2
    unsigned long long temphv1=0;
    unsigned long long temphv2=0;

    for (int i=0; i<startNode->Nchildren; ++i) {
    	unsigned long long t1 = temphv1;
    	unsigned long long t2 = temphv2;
    	unsigned long long h1 = startNode->child[i]->hv1;
    	unsigned long long h2 = startNode->child[i]->hv2;

    	if ( add_of(temphv1, t1, h1) ) {
    		cout << "ERROR: ullong add overflow!!!\n";
    		cout << "t1=" << t1 << " h1=" << h1 << " t1+h1=" << t1+h1 << endl;
    		exit(0);
    	}
    	if ( add_of(temphv2, t2, h2) ) {
    		cout << "ERROR: ullong add overflow!!!\n";
    		cout << "t2=" << t2 << " h2=" << h2 << " t2+h2=" << t2+h2 << endl;
    		exit(0);
    	}
	  }

		// Check overflow
		unsigned long long temp1 = temphv1 % m1;
		unsigned long long temp2 = temphv2 % m2;
    	startNode->hv1 = temp1;
    	startNode->hv2 = temp2;

    // Store bitstrings in hash table
    if (numBitstr < NUM_TAXA-2) {
      vvec_hashrf.hashing_bs_without_type2_nbits(treeIdx, NUM_TAXA, startNode->hv1, startNode->hv2, dist, WEIGHTED);   // without TYPE-III using n-bits (Hash-RF)
    }
  }
}



static void
print_rf_short_matrix(
  vector< vector<unsigned short> > &SHORTSIM,
  unsigned options,
  string outfile)
{
  ofstream fout;
  if (outfile != "")
  {
    fout.open(outfile.c_str());
  }

  switch (options) {
  case 0:
    return;
  case 1:
    cout << "\nRobinson-Foulds distance (list format):\n";

    if (outfile == "") {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          if (i == j)
            cout << "<" << i << "," << j << "> " << 0 << endl;
          else {
              cout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << endl;
          }
        }
      }
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          if (i == j)
            fout << "<" << i << "," << j << "> " << 0 << endl;
          else {
              fout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << endl;
          }
        }
      }
    }

    break;
  case 2:
    cout << "\nRobinson-Foulds distance (rate):\n";
    if (WEIGHTED) {cout << "Fatal error: RF rate is only for unweighted RF distance.\n"; exit(0);}

    if (outfile == "") {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          cout << "<" << i << "," << j << "> ";
          if (i==j)
            cout << 0 << endl;
          else
            cout << (float) ((NUM_TAXA-3)-((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i) {
        for (unsigned int j = 0; j < NUM_TREES; ++j) {
          fout << "<" << i << "," << j << "> ";
          if (i==j)
            fout << 0 << endl;
          else
            fout << (float) ((NUM_TAXA-3)-((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
        }
      }
    }
    break;
  case 3:
    cout << "\nRobinson-Foulds distance (matrix format):\n";
    if (outfile == "") {
      for (unsigned int i = 0; i < NUM_TREES; ++i)  {
        for (unsigned int j = 0; j < NUM_TREES; ++j)  {
//	        for (unsigned int j = i; j < NUM_TREES; ++j)  {
          if (i == j)
            cout << "0" << ' ';
          else
            cout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << ' ';
        }
        cout << endl;
      }
      cout << endl;
    }
    else {
      for (unsigned int i = 0; i < NUM_TREES; ++i)  {
        for (unsigned int j = 0; j < NUM_TREES; ++j)  {
          if (i == j)
            fout << "0" << ' ';
          else
            fout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << ' ';
        }
        fout << endl;
      }
      fout << endl;
    }
    break;
  case 4:
  	if (outfile == "") {
	    for (size_t i = 0; i < NUM_TREES; ++i) {
		    for	(size_t j = 0; j < i; ++j)
		      std::cout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << " ";
		    std::cout << "\n";
		 	}
		}
		else {
			for (size_t i = 0; i < NUM_TREES; ++i) {
		    for	(size_t j = 0; j < i; ++j)
		      fout << (NUM_TAXA-3)-(float((SHORTSIM[i][j] + SHORTSIM[j][i])/2)) << " ";
		    	fout << "\n";
		 	}
		}
    break;
  }

  if (outfile != "")
    fout.close();
}

static void
print_rf_float_matrix(
	vector< vector<float> > &SIM,
	unsigned options,
	string outfile)
{
	ofstream fout;
	if (outfile != "") {
		fout.open(outfile.c_str());
	}

	switch (options) {
		case 0:
			return;
		case 1:
			cout << "\nRobinson-Foulds distance (list format):\n";

				if (outfile == "") {
					for (unsigned i = 0; i < NUM_TREES; ++i) {
						for (unsigned j = 0; j < NUM_TREES; ++j) {
							if (i == j)
								cout << "<" << i << "," << j << "> " << 0 << endl;
							else {
								if (!WEIGHTED)
									cout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << endl;
								else
									cout << "<" << i << "," << j << "> " << float((SIM[i][j] + SIM[j][i])/4) << endl;
							}
						}
					}
				}
				else {
					for (unsigned i = 0; i < NUM_TREES; ++i) {
						for (unsigned j = 0; j < NUM_TREES; ++j) {
							if (i == j)
								fout << "<" << i << "," << j << "> " << 0 << endl;
							else {
								if (!WEIGHTED)
									fout << "<" << i << "," << j << "> " << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << endl;
								else
									fout << "<" << i << "," << j << "> " << float((SIM[i][j] + SIM[j][i])/4) << endl;
							}
						}
					}
				}

			break;
		case 2:
			cout << "\nRobinson-Foulds distance (rate):\n";
			if (WEIGHTED) {
					cout << "Fatal error: RF rate is only for unweighted RF distance.\n";
					exit(0);
			}

			if (outfile == "") {
				for (unsigned i = 0; i < NUM_TREES; ++i) {
					for (unsigned j = 0; j < NUM_TREES; ++j) {
						cout << "<" << i << "," << j << "> ";
						if (i==j)
							cout << 0 << endl;
						else
							cout << (float) ((NUM_TAXA-3)-((SIM[i][j] + SIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
					}
				}
			}
			else {
				for (unsigned i = 0; i < NUM_TREES; ++i) {
					for (unsigned j = 0; j < NUM_TREES; ++j) {
						fout << "<" << i << "," << j << "> ";
						if (i==j)
							fout << 0 << endl;
						else
							fout << (float) ((NUM_TAXA-3)-((SIM[i][j] + SIM[j][i])/2)) / (NUM_TAXA-3) * 100 << endl;
					}
				}
			}
			break;
		case 3:
			cout << "\nRobinson-Foulds distance (matrix format):\n";
			if (outfile == "") {
				for (unsigned i = 0; i < NUM_TREES; ++i)  {
					for (unsigned j = 0; j < NUM_TREES; ++j)  {
//	        for (unsigned j = i; j < NUM_TREES; ++j)  {
						if (i == j)
							cout << " 0 " << ' ';
						else
							if (WEIGHTED)
								cout << float((SIM[i][j] + SIM[j][i])/4) << ' ';
							else
								cout << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << ' ';
					}
					cout << endl;
				}
				cout << endl;
			}
			else {
				for (unsigned i = 0; i < NUM_TREES; ++i)  {
					for (unsigned j = 0; j < NUM_TREES; ++j)  {
						if (i == j)
							fout << " 0 " << ' ';
						else
							if (WEIGHTED)
								fout << float((SIM[i][j] + SIM[j][i])/4) << ' ';
							else
								fout << (NUM_TAXA-3)-(float((SIM[i][j] + SIM[j][i])/2)) << ' ';
					}
					fout << endl;
				}
				fout << endl;
			}
			break;
	}

	if (outfile != "")
		fout.close();
}

int main(int argc, char** argv)
{
  string outfilename;
  string infilename;
  bool bUbid = false; // for counting the number of unique bipartitions

  // TCLAP
  try {

    // Define the command line object.
    string 	helpMsg  = "HashRF\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of HashRF only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";

    helpMsg += "Print out mode: (defualt = matrix).\n";
    helpMsg += "   -p no, no output (default).\n";
    helpMsg += "   -p list, print RF distance in list format.\n";
    helpMsg += "   -p rate, print RF distance rate in list format.\n";
    helpMsg += "   -p matrix, print reuslting distance in matrix format.\n";

    helpMsg += "File option: \n";
    helpMsg += "   -o <export-file-name>, save the distance result in a file.\n";

    helpMsg += "Weighted RF distance: to select RF distance mode between weighted and unweighted (defualt = unweighted).\n";
    helpMsg += "   -w, compute weighted RF distance.\n";

    helpMsg += "Specify c value: \n";
    helpMsg += "   -c <rate>, specify c value (default: 1000) \n";

    helpMsg += "Examples: \n";
    helpMsg += "  hashf foo.tre 1000\n";
    helpMsg += "  hashf bar.tre 1000 -w\n";
    helpMsg += "  hashf bar.tre 1000 -w -p matrix\n";
    helpMsg += "  hashf bar.tre 1000 -w -p list\n";
    helpMsg += "  hashf bar.tre 1000 -w -p list -o output.dat\n";

    TCLAP::CmdLine cmd(helpMsg, ' ', "6.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", true, 2, "Number of trees"  );
    cmd.add( numtreeArg );

    TCLAP::SwitchArg weightedSwitch("w", "weighted", "Compute weighted RF distance", false);
    cmd.add( weightedSwitch );

    TCLAP::ValueArg<string> printArg("p", "printoptions", "print options", false, "matrix", "Print options");
    cmd.add( printArg );

    TCLAP::ValueArg<unsigned int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    cmd.add( cArg );

    TCLAP::SwitchArg ubidSwitch("u", "uBID", "unique BID count", false);
    cmd.add( ubidSwitch );

    TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "", "Output file name");
    cmd.add( outfileArg );

    cmd.parse( argc, argv );

    NUM_TREES = numtreeArg.getValue();

    if (NUM_TREES == 0) {
      string strFileLine;
      unsigned long ulLineCount;
      ulLineCount = 0;

      ifstream infile(argv[1]);
      if (infile) {
        while (getline(infile, strFileLine)) {
          ulLineCount++;
        }
      }
      cout << "*** Number of trees in the input file: " << ulLineCount << endl;
      NUM_TREES = ulLineCount;

      infile.close();
    }

    if (NUM_TREES < 2) {cerr << "Fatal error: at least two trees expected.\n"; exit(2);}

    if (weightedSwitch.getValue())
      WEIGHTED = true;

    if (printArg.getValue() != "matrix") {
      string printOption = printArg.getValue();
      if (printOption == "no")
        PRINT_OPTIONS = 0;
      if (printOption == "list")
        PRINT_OPTIONS = 1;
      if (printOption == "rate")
        PRINT_OPTIONS = 2;
      if (printOption == "cmb")
        PRINT_OPTIONS = 4;
    }

    if (cArg.getValue())
      C = cArg.getValue();

    if (ubidSwitch.getValue())
      bUbid = ubidSwitch.getValue();

    outfilename = outfileArg.getValue();

  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }


  /*****************************************************/
//  cout << "*** Reading a tree file and parsing the tree for taxon label collection ***\n";
  /*****************************************************/

  NEWICKTREE *newickTree;
  int err;
  FILE *fp;
  fp = fopen(argv[1], "r");
  if (!fp) { cout << "Fatal error: file open error\n"; exit(0); }

  newickTree = loadnewicktree2(fp, &err);
  if (!newickTree) {
    switch (err) {
    case -1:
      printf("Out of memory\n");
      break;
    case -2:
      printf("parse error\n");
      break;
    case -3:
      printf("Can't load file\n");
      break;
    default:
      printf("Error %d\n", err);
    }
  }

  /*****************************************************/
  cout << "\n*** Collecting the taxon labels ***\n";
  /*****************************************************/
  LabelMap lm;

  try	{
    GetTaxaLabels(newickTree->root, lm);
  }
  catch (LabelMap::AlreadyPushedEx ex) {
    cerr << "Fatal error: The label '" << ex.label << "' appeard twice in " << endl;
    exit(2);
  }
  NUM_TAXA = lm.size();
  cout << "    Number of taxa = " << lm.size() << endl;
  killnewicktree(newickTree);
  fclose(fp);


  /*****************************************************/
  cout << "\n*** Reading tree file and collecting bipartitions ***\n";
  /*****************************************************/
  HashRFMap vvec_hashrf; // Class HashRFMap

  ////////////////////////////
  // Init hash function class
  ////////////////////////////
  unsigned long long M1=0;
  unsigned long long M2=0;

 	vvec_hashrf.uhashfunc_init(NUM_TREES, NUM_TAXA, C);

  M1 = vvec_hashrf._HF.getM1();
  M2 = vvec_hashrf._HF.getM2();
  vvec_hashrf._hashtab2.resize(M1);

  fp = fopen(argv[1], "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}

  for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {

    newickTree = loadnewicktree2(fp, &err);
    if (!newickTree) {
      switch (err) {
      case -1:
        printf("Out of memory\n");
        break;
      case -2:
        printf("parse error\n");
        break;
      case -3:
        printf("Can't load file\n");
        break;
      default:
        printf("Error %d\n", err);
      }
    }
    else {
      unsigned int numBitstr=0;

      dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2);

      killnewicktree(newickTree);
    }
  }

  cout << "    Number of trees = " << NUM_TREES << endl;
  fclose(fp);


  /*****************************************************/
  cout << "\n*** Compute distance ***\n";
  /*****************************************************/
//  vector< vector<unsigned short> > SIM(NUM_TREES, vector<unsigned short>(NUM_TREES, 0)); // similarity matrix
  typedef vector< vector<float> > FLOAT_MATRIX_T;
  typedef vector< vector<unsigned short> > SHORT_MATRIX_T;
  SHORT_MATRIX_T SHORTSIM;
  FLOAT_MATRIX_T FLOATSIM;

  //---------------------------------------- UNWEIGHTED --------------------------------------------
  if (!WEIGHTED) //  for unweighted RF
    SHORTSIM = SHORT_MATRIX_T (NUM_TREES, vector<unsigned short>(NUM_TREES,0));
  else // for weighted RF
    FLOATSIM = FLOAT_MATRIX_T (NUM_TREES, vector<float>(NUM_TREES,0.0));

	unsigned long uBID = 0;

	if (!WEIGHTED) {    // unweighted
    for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
      unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size();

			if (sizeVec) {

      	uBID += sizeVec;

      	if (!bUbid) {
  	      for (unsigned int i=0; i<sizeVec; ++i) {
  	        unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();

  	        if (sizeTreeIdx > 1) {
  	          for (unsigned int j=0; j<sizeTreeIdx; ++j) {
  	            for (unsigned int k=0; k<sizeTreeIdx; ++k) {
  	              if (j == k) continue;
  	              else {
  	                SHORTSIM[vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j]][vvec_hashrf._hashtab2[hti][i]._vec_treeidx[k]] += 1;
  	              }
  	            }
  	          }
  	        }
  	      }
	      } //if
	    }
    }
  }
  //---------------------------------------- WEIGHTED --------------------------------------------
  else {
    vvec_hashrf._hashtab.resize(vvec_hashrf._hashtab2.size());
    for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
      unsigned int sizeLinkedList = vvec_hashrf._hashtab2[hti].size();
      if (sizeLinkedList > 0) {
        for (unsigned int i1=0; i1<sizeLinkedList; ++i1) {
          unsigned int bidi = vvec_hashrf._hashtab2[hti][i1]._vec_treeidx.size();
          for (unsigned int i2=0; i2<bidi; ++i2) {
            BUCKET_STRUCT_T bk;
            bk._hv2 = vvec_hashrf._hashtab2[hti][i1]._hv2;
            bk._t_i = vvec_hashrf._hashtab2[hti][i1]._vec_treeidx[i2];
            bk._dist = vvec_hashrf._hashtab2[hti][i1]._vec_dist[i2];
            vvec_hashrf._hashtab[hti].push_back(bk);
          }
        }
      }
    }
    vvec_hashrf._hashtab2.clear();

    for (unsigned int hti=0; hti<vvec_hashrf._hashtab.size(); ++hti) {

      unsigned int sizeLinkedList = vvec_hashrf._hashtab[hti].size();

      if (sizeLinkedList > 1) {
        vector<unsigned long> vec_hv2;
        vector<unsigned long>::iterator itr_vec_hv2;

        ////////////////////////////////////////////////
        // Collect unique hv2 values in the linked list
        ////////////////////////////////////////////////
        for (unsigned int i=0; i<sizeLinkedList; ++i) {
          unsigned long hv2 = vvec_hashrf._hashtab[hti][i]._hv2;
          if (vec_hv2.empty())
            vec_hv2.push_back(hv2);
          else {
            itr_vec_hv2 = find(vec_hv2.begin(), vec_hv2.end(), hv2);
            if (itr_vec_hv2 == vec_hv2.end())
              vec_hv2.push_back(hv2);
          }
        }

        // distance array
        vector< vector<float> > vvec_dist(vec_hv2.size(), vector<float>(NUM_TREES, 0));

        /////////////////////////////////////////////////////////////
        // SET THE distance array with distance at proper tree index
        /////////////////////////////////////////////////////////////
        for (unsigned int i=0; i<sizeLinkedList; ++i) {
          for (unsigned int j=0; j<vec_hv2.size(); ++j) {
            if (vvec_hashrf._hashtab[hti][i]._hv2 == vec_hv2[j])
              vvec_dist[j][vvec_hashrf._hashtab[hti][i]._t_i] = vvec_hashrf._hashtab[hti][i]._dist;
          }
        }

        /////////////////////////////////////
        // UPDATE FLOATSIM MATIRX USING vvec_dist
        /////////////////////////////////////
        for (unsigned int i=0; i<vvec_dist.size(); ++i) {

          for (unsigned int j=0; j<vvec_dist[i].size(); ++j) {
            for (unsigned int k=0; k<vvec_dist[i].size(); ++k) {
              if (j == k) continue;
              else
                FLOATSIM[j][k] += abs(vvec_dist[i][j] - vvec_dist[i][k]);
            }
          }
        }

        vec_hv2.clear();
        vvec_dist.clear();
      }
      else if (sizeLinkedList == 1) {

        /////////////////////////////////////////////////////
        // PROPAGATE the dist value TO OTHER TREES' distance
        /////////////////////////////////////////////////////
        for (unsigned int i=0; i<NUM_TREES; ++i) {

          if (i == vvec_hashrf._hashtab[hti][0]._t_i)
            continue;
          else {
            FLOATSIM[i][vvec_hashrf._hashtab[hti][0]._t_i] += vvec_hashrf._hashtab[hti][0]._dist;
            FLOATSIM[vvec_hashrf._hashtab[hti][0]._t_i][i] += vvec_hashrf._hashtab[hti][0]._dist;
          }
        }
      }
    }
  }

  cout << "    # of unique BIDs = " << uBID << endl;

	if (!WEIGHTED)
  	print_rf_short_matrix(SHORTSIM, PRINT_OPTIONS, outfilename);
  else
  	print_rf_float_matrix(FLOATSIM, PRINT_OPTIONS, outfilename);


  /*****************************************************/
//  cout << "\n*** Print statistics ***\n";
  /*****************************************************/
  // CPU time comsumed
  struct rusage a;
  if (getrusage(RUSAGE_SELF,&a) == -1) { cerr << "Fatal error: getrusage failed.\n";  exit(2); }
  cout << "\n    Total CPU time: " << a.ru_utime.tv_sec+a.ru_stime.tv_sec << " sec and ";
  cout << a.ru_utime.tv_usec+a.ru_stime.tv_usec << " usec.\n";


  /*****************************************************/
//  cout << "\n*** Clear allocated memory***\n";
  /*****************************************************/
  SHORTSIM.clear();
  FLOATSIM.clear();
  vvec_hashrf.hashrfmap_clear();


  return 1;
}


// eof
