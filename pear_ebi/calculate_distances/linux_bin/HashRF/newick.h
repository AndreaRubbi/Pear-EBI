
#ifndef newick_h
#define newick_h

typedef struct newicknode
{
  int Nchildren;              /* number of children (0 for leaves) */
  char *label;                /* node label, can be null */
  double weight;              /* node weight */
  struct newicknode **child;  /* list of children */
  unsigned long long hv1;
  unsigned long long hv2;
} NEWICKNODE;

typedef struct
{
  NEWICKNODE *root;
} NEWICKTREE;

NEWICKTREE *loadnewicktree(char *fname, int *error);
NEWICKTREE *loadnewicktree2(FILE *fp, int *error);
NEWICKTREE *floadnewicktree(FILE *fp, int *error);
void killnewicktree(NEWICKTREE *tree);
char *makenewicklabel(const char *str);
void printnewicktree(NEWICKTREE *tree);

#endif
