/*
  Newick format tree loader
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

static void killnoder(NEWICKNODE *node);
static NEWICKNODE *loadnode(FILE *fp, int *error);
static int addchild(NEWICKNODE *parent, NEWICKNODE *child);
static NEWICKNODE *loadleaf(FILE *fp, int *error);
static char *readlabelandweight(FILE *fp, double *weight, int *error);
static char *readlabel(FILE *fp);
static char *readquotedlabel(FILE *fp);
static char *readplainlabel(FILE *fp);
static void skipspace(FILE *fp);
static char *mystrdup(const char *str);
static int mystrcount(const char *str, int ch);

/*
  load a  newick tree from a file
  Params; fname - path to file
          error - return for error code
  Returns: the loaed object on success, 0 on fail
  Error codes: -1 out of memory
               -2 parse error
               -3 can't open file
 */
NEWICKTREE *loadnewicktree(char *fname, int *error)
{
  NEWICKTREE *answer;
  FILE *fp;
  int err;

  fp = fopen(fname, "r");
  if(!fp)
  {
    err = -3;
    goto error_exit;
  }
  answer = floadnewicktree(fp, &err);
  if(!answer)
    goto error_exit;
  if(error)
    *error = 0;
  return answer;

 error_exit:
  if(error)
    *error = err;
  return 0;

}


NEWICKTREE *loadnewicktree2(FILE *fp, int *error)
{
  NEWICKTREE *answer;
//  FILE *fp;
  int err;

//  fp = fopen(fname, "r");
//  if(!fp)
//  {
//    err = -3;
//    goto error_exit;
//  }
  answer = floadnewicktree(fp, &err);
  if(!answer)
    goto error_exit;
  if(error)
    *error = 0;
  return answer;

 error_exit:
  if(error)
    *error = err;
  return 0;

}

/*
  load newick tree from an opened stream
  Params: fp - input stream
          error - return for error
  Returns: pointer to object, 0 on fail
  Error codes -1 out of memory
              -2 parse error
  Notes: format allows several trees to be stroed in a file
 */
NEWICKTREE *floadnewicktree(FILE *fp, int *error)
{
  NEWICKTREE *answer;
  int err;
  int ch;

  answer = malloc(sizeof(NEWICKTREE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch == '(' )
  {
    answer->root = loadnode(fp, &err);
    if(!answer->root || err != 0)
      goto error_exit;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch != ';')
  {
    err = -2;
    goto error_exit;
  }

  if(error)
    *error = 0;
  return answer;

 error_exit:
  killnewicktree(answer);
  if(error)
    *error = err;
  return 0;

}

/*
  newick tree destructor
 */
void killnewicktree(NEWICKTREE *tree)
{
  if(tree)
  {
    killnoder(tree->root);
    free(tree);
  }
}

/*
  turns a string into a label suitable for use in the tree
  Params: str - the string
  Returns: modified string, 0 on out of memory
  Notes: strings containing spaces have them replaced by underscores
         strings contianing illegal characters are quoted
         null pointer is returned as the empty string ""
 */
char *makenewicklabel(const char *str)
{
  const char *cptr;
  char *vptr;
  int needsquote = 0;
  char *answer;

  if(!str)
    return mystrdup("");
  cptr = str;
  while(*cptr)
  {
    if(strchr("\'()[]:;,", *cptr))
      needsquote = 1;
    if(isspace(*cptr) && *cptr != ' ')
      needsquote = 1;
    cptr++;
  }
  if(needsquote)
  {
    answer = malloc(strlen(str) + 2 + mystrcount(str, '\''));
    vptr = answer;
    *vptr++ = '\'';
    while(*str)
    {
      *vptr = *str;
      if(*str == '\'')
      {
        vptr++;
        *vptr = '\'';
      }
      vptr++;
      str++;
    }
    *vptr++ = '\'';
    *vptr = 0;
  }
  else
  {
    answer = mystrdup(str);
    vptr = answer;
    while(*vptr)
    {
      if(*vptr == ' ' )
 	*vptr = '_';
      vptr++;
    }
  }

  return answer;
}

/*
  node destructor (recursive)
 */
static void killnoder(NEWICKNODE *node)
{
  int i;

  if(!node)
    return;

  for(i=0;i<node->Nchildren;i++)
  {
    killnoder(node->child[i]);
  }
  free(node->label);
  free(node->child);
  free(node);
}

/*
  load a node from the file
  Params: fp - the input stream
          error - return for error
  Returns: node loaded.
  Notes: recursive. Expects the opening parenthesis to have been eaten
 */
static NEWICKNODE *loadnode(FILE *fp, int *error)
{
  NEWICKNODE *answer;
  int err;
  NEWICKNODE *child = 0;
  int ch;

  answer = malloc(sizeof(NEWICKNODE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }

  answer->Nchildren = 0;
  answer->label = 0;
  answer->child = 0;
  answer->hv1 = 0;
  answer->hv2 = 0;

  skipspace(fp);
  do
  {
    ch = fgetc(fp);
    if(ch == '(')
    {
      child = loadnode(fp, &err);
      if(!child)
        goto error_exit;

      if( addchild(answer, child ) == -1)
      {
        err = -1;
        goto error_exit;
      }
      child = 0;
    }
    else
    {
      ungetc(ch, fp);
      child = loadleaf(fp, &err);
      if(!child)
        goto error_exit;

      if(addchild(answer, child) == -1)
      {
        err = -1;
        goto error_exit;
      }
      child = 0;
    }
    skipspace(fp);
    ch = fgetc(fp);
  } while(ch == ',');

  if(ch == ')')
  {
    answer->label = readlabelandweight(fp, &answer->weight, &err);
    if(err)
      goto error_exit;
  }
  else
  {
    err = -2;
    goto error_exit;
  }

  if(error)
    *error = 0;
  return answer;
 error_exit:
  if(child)
    killnoder(child);
  killnoder(answer);
  if(error)
    *error = err;
  return 0;
}

/*
  add a child to a node
  Params: parent - parent node
          child - child to add
  Returns: 0 on siccess, -1 on fail
 */
static int addchild(NEWICKNODE *parent, NEWICKNODE *child)
{
  NEWICKNODE **temp;

  temp = realloc(parent->child, (parent->Nchildren + 1) * sizeof(NEWICKNODE *));
  if(!temp)
    return -1;
  parent->child = temp;
  parent->child[parent->Nchildren] = child;
  parent->Nchildren++;

  return 0;
}


/*
  load a leaf node
  Params: fp - the input stream
          error - return for error
  Returns: node object

 */
static NEWICKNODE *loadleaf(FILE *fp, int *error)
{
  NEWICKNODE *answer;
  int err;

  answer = malloc(sizeof(NEWICKNODE));
  if(!answer)
  {
    err = -1;
    goto error_exit;
  }

  answer->Nchildren = 0;
  answer->child = 0;
  answer->label = readlabelandweight(fp, &answer->weight, &err);
  if(err)
    goto error_exit;
  return answer;

  error_exit:
  if(error)
    *error = err;
  if(answer)
  {
    free(answer->label);
    free(answer);
  }
  return 0;
}

/*
  read label, colon and associated weight
  Params: fp - the input stream
          weight - return for weight
          error - return for error
  Returns: the label
  Notes: a null label is not an error
 */
static char *readlabelandweight(FILE *fp, double *weight, int *error)
{
  char *answer;
  int ch;

  *weight = 0;
  answer = readlabel(fp);
  if(!answer)
  {
    *error = -1;
    return 0;
  }
  skipspace(fp);
  ch = fgetc(fp);
  if(ch == ':')
  {
    if( fscanf(fp, "%lf", weight) != 1)
    {
      *error = -2;
      free(answer);
      return 0;
    }
  }
  else
    ungetc(ch, fp);

  if(*answer == 0)
  {
    free(answer);
    answer = 0;
  }
  *error = 0;
  return answer;
}

/*
  readlabel - read a label from stream
  Params: fp - the input stream
  Returns: allocarted label
  Notes: null label is the null string ""
 */
static char *readlabel(FILE *fp)
{
  char *answer;
  int ch;

  skipspace(fp);
  ch = fgetc(fp);
  ungetc(ch, fp);
  if(ch == '\'')
    answer = readquotedlabel(fp);
  else
    answer = readplainlabel(fp);

  return answer;
}

/*
  read a quoted label from stream
  Params: fp - the input stream
  Returns: the label. null string is ""
 */
static char *readquotedlabel(FILE *fp)
{
  char *answer;
  char *temp;
  int capacity = 32;
  int len = 0;
  int ch;

  answer = malloc(capacity);
  if(!answer)
    return 0;
  ch = fgetc(fp);
  while( (ch = fgetc(fp)) != EOF)
  {
    if(ch == '\'')
    {
      ch = fgetc(fp);
      if(ch == '\'')
        answer[len++] = (char) ch;
      else
      {
        ungetc(ch, fp);
        answer[len] = 0;
        return answer;
      }
    }
    else
      answer[len++] = (char) ch;
    if(len == capacity - 1)
    {
      temp = realloc(answer, capacity * 2);
      if(!temp)
      {
        free(answer);
        return 0;
      }
      answer = temp;
      capacity = capacity * 2;
    }
  }
  answer[len] = 0;
  return answer;
}

/*
  read aplian (unquoted) label
  Params: fp - the open file
  Returns: the label, null string is ""
 */
static char *readplainlabel(FILE *fp)
{
  char *answer;
  char *temp;
  int len = 0;
  int capacity = 32;
  int ch;

  answer = malloc(capacity);
  if(!answer)
    return 0;

  while( (ch = fgetc(fp)) != EOF)
  {
    if(isspace(ch) || strchr(" ()[]\':;,", ch))
    {
      ungetc(ch, fp);
      answer[len] = 0;
      return answer;
    }
    if(ch == '_')
      ch = '_'; // ssj
    answer[len++] = (char) ch;
    if(len == capacity - 1)
    {
      temp = realloc(answer, capacity * 2);
      if(!temp)
      {
        free(answer);
        return 0;
      }
      answer = temp;
      capacity = capacity * 2;
    }
  }

  answer[len] = 0;
  return answer;
}

/*
  consume space in input stream.
  Params: fp - the input stream
 */
static void skipspace(FILE *fp)
{
  int ch;

  while( (ch = fgetc(fp)) != EOF )
    if(!isspace(ch))
    {
      ungetc(ch, fp);
      break;
    }
}

/*
  strduo implemetation
 */
static char *mystrdup(const char *str)
{
  char *answer;

  answer = malloc(strlen(str) + 1);
  if(answer)
    strcpy(answer, str);

  return answer;
}

/*
  count the number of instances of ch in the string
  Params: str - serach string
          ch - character.
  Returns: number of times ch appears in str.
 */
static int mystrcount(const char *str, int ch)
{
  int answer = 0;

  while(*str)
    if(*str++ == ch)
      answer++;

  return answer;
}

/*
  print out a node (test function)
  Parmas: node - node to print
          indent - indent level
  Notes: recursive
*/
static void printnewicknode(NEWICKNODE *node, int indent)
{
  int i;

  for(i=0;i<indent;i++)
    printf(" ");
  printf("%s : %f\n", node->label, node->weight);
  for(i=0;i<node->Nchildren;i++)
    printnewicknode(node->child[i], indent + 5);
}

/*
  print the Newick tree in a human-readable format (test function)
  Params: tree - the tree
 */
void printnewicktree(NEWICKTREE *tree)
{
  printnewicknode(tree->root, 0);
}

int newickmain(int argc, char **argv)
{
  NEWICKTREE *tree;
  int err;
  char *label;

  if(argc != 2)
    fprintf(stderr, "Newick tree loader\n");
  else
  {
    tree = loadnewicktree(argv[1], &err);
    if(!tree)
    {
      switch(err)
      {
      case -1: printf("Out of memory\n"); break;
      case -2: printf("parse error\n"); break;
      case -3: printf("Can't load file\n"); break;
      default:
        printf("Error %d\n", err);
      }
    }
    else
    {
      printf("Loaded\n");
      printnewicktree(tree);
    }
    killnewicktree(tree);
  }
  label = makenewicklabel("My name is \'Fred()\'");
  printf("***%s***\n", label);
  free(label);

  return 0;
}
