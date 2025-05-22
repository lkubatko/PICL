#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#define MAXPASS 5
#define SMAXPASS 5
#define DELTA 0.00001
#define SDELTA 0.0001
#define FTOL 0.00001

extern int ntaxa, nspecies, nsite, nquarts, num_unique_quarts, num_unique, include_gaps, num_no_gaps, verbose, model, ncat, random_tree;
extern int anneal, anneal_bl, user_bl, max_it, max_it_bl, seedj, seedk;

extern double ci, max_cl;
extern float theta, beta, mu, ratepar, invpar;

extern int *parents, *parents_temp, *ppTwoRow[2], *ppTwoRow_temp[2], *ppTwoRow_best[2], *ppTwoRowQuart[2], *filled_ind, *seq_counter, *qvec;
extern int **ppBase_full, **ppBase, **ppBase_unique, *site_counter, **ppSp_assign, **ppNodeChildren, **ppNodeChildrenLeftQuart, **ppNodeChildrenRightQuart;
extern int max_it;
extern double *TimeVec, *TimeVec_temp, *TimeVec_init, *TimeVec_best, *rvals, *TimeVecQuart, **ppLengthMat, **ppMatrix;
extern double smat[10][10],amat[12][12];

extern FILE *out;

typedef char naym[12];

struct Node {
  char *iTree;
  int first_hit;
  int first_anneal;
  double max_lik;
  int times_hit;
  int loc_ppTR;
  struct Node *Next;
};

typedef struct Node *Link;


struct CLquart {
  int spprobs[15];
  int ncherries;
  double *t1;
  double *t2;
  double *t3;
};

typedef struct CLquart CLq;

extern Link Head, currenttree;
extern naym *taxname, *speciesname, *psname;
extern CLq **StoreQuarts;

/* in main.c */
//long int random_seed();
int MemAlloc();
void transfer(char **cbase, int **nbase, int pres_loc, int trans_length);
naym* getSequenceData(FILE *infile);
void unique_sites();

/* in randomtree.c */
int find_parentr(int target);
int find_genr(int current_node);
void build_random_tree(int seed1, int seed2);
void assign_times();

/* in tree.c */
int FindParentI(int target);
int NextNode(int next);
double ReadLength(FILE *fp);
int AddToNode(int current,int next,int last_tip,FILE *fp,double read_length);
int FindLastOpen(int current);
int CloseBack(int open, int current, int last_char, int last_tipp, FILE *fp, double read_length);
double ReadTree(FILE *fp);
double FindTotalTime();
void CalcTimeVec(double total, double rate);

/* in compare.c */
void write_species_tree(int node, int previous_node);
void print_PHYLIP();
void remove_CONSTANT();
unsigned long long binomial(int n, int k);
unsigned long long subset_to_index(int a, int b, int c, int d);
void index_to_subset(unsigned long long index, int *a, int *b, int *c, int *d, int n);

/* in complik.c */
int find_parent(int target);
int IsCherry(int q, int r);
int IsCherryQuartet(int q, int r, int w, int x);
int CountCherries(int q, int r, int s, int t);
int find_gen(int current_node);
void make_indmat();
int FindMRCA(int q, int r);
void AllMyChildrenLeftQuartLeft(int mynode);
void AllMyChildrenRightQuartetLeft(int mynode);
void AllMyChildrenLeftQuartRight(int mynode);
void AllMyChildrenRightQuartetRight(int mynode);
void AllMyChildrenQuartet(int mynode);
void FindDupQuarts(int is_symm, int spcounts[15]);
int GetQuartetTree(int tax1, int tax2, int tax3, int tax4, int OrderedVec[5]);
void CountQuartetSitePatterns(int tax1, int tax2, int tax3, int tax4, int pvec[15]);
double SymmetricQuartetLikelihood(int nn);
double AsymmetricQuartetLikelihood(int nn);
double SymmetricQuartetLikelihood_msnp(int nn);
double AsymmetricQuartetLikelihood_msnp(int nn);
double SymmetricQuartetLikelihood_ratevar(int nn);
double AsymmetricQuartetLikelihood_ratevar(int nn);
double GetCompLik();
double GetCompLik_msnp();
double GetCompLik_ratevar();
void GetRateParams();
void swap(int *a, int *b);
void sort4(int arr[]);

/* in nodeopt.c */
void AllMyChildrenLeft(int mynode);
void AllMyChildrenRight(int mynode);
void AllMyChildren(int mynode);
double GetDerCompLik(int mynode);
double GetCompMOM(int ltax1, int ltax2, int ltax3, int ltax4, int issymmetric, int whichtau);
double GetMomentEstimate(int mynode);

/* in anneal.c */
void anneal_full();
void bl_anneal_full();
void bl_uphill_full();
void bl_anneal_msnp(); 
void bl_uphill_msnp();
void bl_anneal_ratevar();
void bl_uphill_ratevar();

/* in boot.c */
void boot_times(int nrep);

/* in opt.c */
void OptAlpha();

/* in trbldg.c */
void trbldg();

/* in trbldg_ratevar.c */
void trbldg_ratevar();

/* in trbldg_msnp.c */
void trbldg_msnp();

