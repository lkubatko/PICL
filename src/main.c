#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include "main.h"
#include "linpack.c"
#include "ranlib.h"
#include "ranlib.c"
#include "asa091.c"
#include "com.c"
#include "randomtree.c"
#include "tree.c"
#include "compare.c"
#include "complik.c"
#include "complik_genetree.c"
#include "nodeopt.c"
#include "anneal.c"
#include "opt.c"
#include "trbldg.c"
#include "trbldg_ratevar.c"
#include "trbldg_msnp.c"
#include "trbldg_genetree.c"
#include "complik_popvar.c"
#include "treebldg_popvar.c"
#include "anneal_popvar.c"
#include "boot.c"
#include "boot_times.c"


// ntaxa = number of species
// nseq = number of individuals
int ntaxa, nspecies, nseq, nsite, nquarts, num_unique_quarts, num_unique, include_gaps, num_no_gaps, verbose, model, ncat, random_tree;
int anneal, anneal_bl, user_bl, max_it, mult_iter, num_reject, max_it_bl, test_increment, seedj, seedk;
int *parents, *parents_temp, *ppTwoRow[2], *ppTwoRow_temp[2], *ppTwoRow_best[2], *ppTwoRow_best_orig[2], *ppTwoRowQuart[2], *filled_ind, *seq_counter, *qvec, *site_counter;
int **ppBase_full, **ppBase, **ppBase_unique, **ppSp_assign, **ppNodeChildren, **ppNodeChildrenLeftQuart, **ppNodeChildrenRightQuart;
int pattern_index[16][16];
double ci, max_cl, curr_anneal_lik, b1opt, prob_bound;
float alpha, lambda, theta, beta, mu, ratepar, invpar;
double *TimeVec, *TimeVec_temp, *TimeVec_init, *TimeVec_best, *TimeVec_best_orig, *TimeVecQuart, *rvals, **ppLengthMat, **ppMatrix;
double smat[10][10],amat[12][12];
double base_weight_table[15][4]; 
FILE *out, *pt, *boot;

Link Head, currenttree;
naym *taxname=NULL, *speciesname=NULL, *psname=NULL;
naym **assignment_vector=NULL;
CLq** StoreQuarts;

int MemAlloc() {

  int i, j;

  TimeVec = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVec==NULL)
    {
      printf("Can't memalloc TimeVec\n");
      return 0;
    }

  TimeVec_temp = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVec_temp==NULL)
    {
      printf("Can't memalloc TimeVec_temp\n");
      return 0;
    }

  TimeVec_init = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVec_init==NULL)
    {
      printf("Can't memalloc TimeVec_init\n");
      return 0;
    }

  TimeVec_best = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVec_best==NULL)
    {
      printf("Can't memalloc TimeVec_best\n");
      return 0;
    }

  TimeVec_best_orig = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVec_best_orig==NULL)
    {
      printf("Can't memalloc TimeVec_best_orig\n");
      return 0;
    }

  TimeVecQuart = (double*)malloc((2*ntaxa+1)*sizeof(double));
  if (TimeVecQuart==NULL)
    {
      printf("Can't memalloc TimeVecQuart\n");
      return 0;
    }  

  rvals = (double*)malloc((ncat+1)*sizeof(double));
  if (rvals==NULL)
    {
      printf("Can't memalloc rvals\n");
      return 0;
    }

  for (i=0; i<2; i++)
    {
      ppTwoRow[i] = (int*)malloc((2*ntaxa)*sizeof(int));
      if (ppTwoRow[i]==NULL)
        {
          printf("Can't memalloc ppTwoRow[%d]\n",i);
          return 0;
        }
    }

  for (i=0; i<2; i++)
    {
      ppTwoRow_temp[i] = (int*)malloc((2*ntaxa)*sizeof(int));
      if (ppTwoRow_temp[i]==NULL)
        {
          printf("Can't memalloc ppTwoRow_temp[%d]\n",i);
          return 0;
        }
    }


  for (i=0; i<2; i++)
    {
      ppTwoRow_best[i] = (int*)malloc((2*ntaxa)*sizeof(int));
      if (ppTwoRow_best[i]==NULL)
        {
          printf("Can't memalloc ppTwoRow_best[%d]\n",i);
          return 0;
        }
    }

  for (i=0; i<2; i++)
    {
      ppTwoRow_best_orig[i] = (int*)malloc((2*ntaxa)*sizeof(int));
      if (ppTwoRow_best_orig[i]==NULL)
        {
          printf("Can't memalloc ppTwoRow_best_orig[%d]\n",i);
          return 0;
        }
    }

  for (i=0; i<2; i++)
    {
      ppTwoRowQuart[i] = (int*)malloc((2*ntaxa)*sizeof(int));
      if (ppTwoRowQuart[i]==NULL)
        {
          printf("Can't memalloc ppTwoRowQuart[%d]\n",i);
          return 0;
        }
    }

  parents = (int*)malloc((2*ntaxa+1)*sizeof(int));
  if (parents==NULL) 
    {
      printf("Can't memalloc parents\n");
      return 0;
    }

  parents_temp = (int*)malloc((2*ntaxa+1)*sizeof(int));
  if (parents_temp==NULL)
   {
      printf("Can't memalloc parents_temp\n");
      return 0;
    }

  filled_ind = (int*)malloc((2*ntaxa)*sizeof(int));

  /* Every use of qvec throughout the codebase resets/reads it with loops
   * shaped "for (i=0; i<nquarts+1; i++)" (main.c, complik.c, anneal.c,
   * boot.c, etc.) -- i.e. it is consistently treated as having nquarts+1
   * valid indices (0..nquarts). This allocation was a genuine,
   * pre-existing off-by-one: nquarts ints instead of nquarts+1, found via
   * AddressSanitizer (see `make debug`) as a heap-buffer-overflow. */
  qvec = (int*)malloc((nquarts+1)*sizeof(int));

  seq_counter = (int*)malloc((ntaxa+1)*sizeof(int));
  if (seq_counter==NULL)
    { 
      printf("Can't memalloc seq_counter\n");
      return 0;
    }	

   ppLengthMat = (double**)malloc((2*ntaxa+1)*sizeof(double*));
  if (ppLengthMat==NULL)
    {
      printf("Can't memalloc ppLengthMat\n");
      return 0;
    }
  
  for (i=0; i<2*ntaxa+1; i++)
    {
      ppLengthMat[i] = (double*)malloc((2*ntaxa+1)*sizeof(double));
      if (ppLengthMat[i]==NULL)
        {
          printf("Can't memalloc ppLengthMat[%d]\n",i);
          return 0;
        }
    }

  ppMatrix = (double**)malloc((nseq+1)*sizeof(double*));
  if (ppMatrix==NULL)
    {
      printf("Can't memalloc ppMatrix\n");
      return 0;
    }

  for (i=0; i<nseq+1; i++)
    {
      /* make_indmat() (complik.c) writes ppMatrix[q][r] for q,r up to and
       * including ntaxa (its loops run q,r < ntaxa+1, i.e. q,r <= ntaxa),
       * so this needs ntaxa+1 doubles (valid indices 0..ntaxa), not ntaxa.
       * The previous ntaxa-sized allocation was a genuine, pre-existing
       * off-by-one heap overflow: it silently corrupted adjacent heap
       * memory without necessarily crashing under glibc's allocator, but
       * reliably segfaults right after tree-reading under stricter
       * allocators (e.g. macOS). Found via AddressSanitizer (see `make
       * debug` in the Makefile). */
      ppMatrix[i] = (double*)malloc((ntaxa+1)*sizeof(double));
      if (ppMatrix[i]==NULL)
        {
          printf("Can't memalloc ppMatrix[%d]\n",i);
          return 0;
        }
    }
  
  
  if(nseq <= 0 || nseq >100000){
    printf("Number of sequences too big or negative");
    exit(1);
  }
  else{
  ppBase_full = (int**)malloc(nseq*sizeof(int*));
  }
  if (ppBase_full==NULL)
    {
      printf("     Can't memalloc ppBase_full.\n");
      return 0;
    }

  for (i=0; i<nseq; i++) {
    ppBase_full[i]=(int*)malloc((nsite+1)*sizeof(int));
    if (ppBase_full[i]==NULL) 
      {
	printf("     Can't memalloc ppBase_full[%d].\n",i);
	return 0;
      }
  }

  if(nseq <= 0 || nseq >100000){
    printf("Number of sequences too big or negative");
    exit(1);
  }
  else{
    ppBase = (int**)malloc(nseq*sizeof(int*));
  }
  
  if (ppBase==NULL)
    {
      printf("     Can't memalloc ppBase.\n");
      return 0;
    }
  
  for (i=0; i<nseq; i++)
    {
      ppBase[i]=(int*)malloc((nsite+1)*sizeof(int));
      if (ppBase[i]==NULL)
	{
	  printf("     Can't memalloc ppBase[%d].\n",i);
	  return 0;
	}
    }
 
  ppBase_unique = (int**)malloc(nseq*sizeof(int*));
  if (ppBase_unique==NULL)
    {
      printf("     Can't memalloc ppBase_unique.\n");
      return 0;
    }
  
  for (i=0; i<nseq; i++)
    {
      ppBase_unique[i]=(int*)malloc((nsite+5)*sizeof(int));
      if (ppBase_unique[i]==NULL)
	{
	  printf("     Can't memalloc ppBase_unique[%d].\n",i);
	  return 0;
	}
    }

  site_counter = (int*)malloc((nsite+1)*sizeof(int));
  if (site_counter==NULL)
    {
      printf("     Can't memalloc site_counter.\n");
      return 0;
    }

  ppSp_assign = (int**)malloc(ntaxa*sizeof(int*));
  if (ppSp_assign==NULL)
    {
      printf("     Can't memalloc ppSp_assign.\n");
      return 0;
    }
  for (i=0; i<ntaxa; i++) 
    { 
      ppSp_assign[i] = (int*)malloc(100*sizeof(int));
      if (ppSp_assign[i]==NULL)
	{
	  printf("     Can't memalloc ppSp_assign[%d].\n",i);
	  return 0;
	}
  }

  ppNodeChildren = (int**)malloc(3*sizeof(int*));
  if (ppNodeChildren==NULL)
   {
     printf("     Can't memalloc ppNodeChildren.\n");
     return 0;
  }
  /* ppNodeChildren, ppNodeChildrenLeftQuart, and ppNodeChildrenRightQuart
   * each store a count in slot [0] followed by up to ntaxa child indices in
   * slots [1..ntaxa] (see the AllMyChildren-family helpers and FindDupQuarts
   * in complik.c, and nodeopt.c, which in the worst case -- every taxon
   * descended from the relevant node -- write index ntaxa). That needs
   * ntaxa+1 ints, not ntaxa; this was a genuine, pre-existing off-by-one
   * heap overflow found via AddressSanitizer (see the debug target in the
   * Makefile). */
  for (i=0; i<3; i++)
   {
     ppNodeChildren[i] = (int*)malloc((ntaxa+1)*sizeof(int));
     if (ppNodeChildren[i]==NULL)
       {
         printf("     Can't memalloc ppNodeChildren[%d].\n",i);
         return 0;
       }
   }

  ppNodeChildrenLeftQuart = (int**)malloc(3*sizeof(int*));
  if (ppNodeChildrenLeftQuart==NULL)
   {
     printf("     Can't memalloc ppNodeChildrenLeftQuart.\n");
     return 0;
  }
  for (i=0; i<3; i++)
   {
     ppNodeChildrenLeftQuart[i] = (int*)malloc((ntaxa+1)*sizeof(int));
     if (ppNodeChildrenLeftQuart[i]==NULL)
       {
         printf("     Can't memalloc ppNodeChildrenLeftQuart[%d].\n",i);
         return 0;
       }
   }

  ppNodeChildrenRightQuart = (int**)malloc(3*sizeof(int*));
  if (ppNodeChildrenRightQuart==NULL)
   {
     printf("     Can't memalloc ppNodeChildrenRightQuart.\n");
     return 0;
  }
  for (i=0; i<3; i++)
   {
     ppNodeChildrenRightQuart[i] = (int*)malloc((ntaxa+1)*sizeof(int));
     if (ppNodeChildrenRightQuart[i]==NULL)
       {
         printf("     Can't memalloc ppNodeChildrenRightQuart[%d].\n",i);
         return 0;
       }
   }

  StoreQuarts = (CLq**)malloc(nquarts*sizeof(CLq*));
  for (i=0; i<nquarts; i++) StoreQuarts[i] = (CLq*)malloc(sizeof(CLq));


  return 1;

}

// Transform character sequences (A,G,C,T) into numerical values (0,1,2,3) 

void transfer(char **cbase, int **nbase, int pres_loc, int trans_length) {

  int i, j;

  for (i=0; i<nseq; i++)  {

    for (j=0; j<trans_length; j++) {

      switch (cbase[i][j]) {
	    
      case 'A': case 'a': nbase[i][j+pres_loc] = 0;  break;
      case 'G': case 'g': nbase[i][j+pres_loc] = 1;  break;
      case 'C': case 'c': nbase[i][j+pres_loc] = 2;  break;
      case 'T': case 't': case 'U': case 'u': nbase[i][j+pres_loc] = 3;  break;
      case 'M': case 'm': nbase[i][j+pres_loc] = 5;  break;
      case 'R': case 'r': nbase[i][j+pres_loc] = 6;  break;
      case 'W': case 'w': nbase[i][j+pres_loc] = 7;  break;
      case 'S': case 's': nbase[i][j+pres_loc] = 8;  break;
      case 'Y': case 'y': nbase[i][j+pres_loc] = 9;  break;
      case 'K': case 'k': nbase[i][j+pres_loc] = 10;  break;
      case 'B': case 'b': nbase[i][j+pres_loc] = 11;  break;
      case 'D': case 'd': nbase[i][j+pres_loc] = 12;  break;
      case 'H': case 'h': nbase[i][j+pres_loc] = 13;  break;
      case 'V': case 'v': nbase[i][j+pres_loc] = 14;  break;
      case 'N': case 'n': nbase[i][j+pres_loc] = 15;  break;
      case '.':nbase[i][j+pres_loc] = nbase[0][j+pres_loc]; break;
      case '-':nbase[i][j+pres_loc]=4; break;
      case '?':nbase[i][j+pres_loc]=15; break;
      default: break;

      }
	  
    }
  }
}



// With number of sequences (nseq), number of sites (nsite),    
// get  DNA sequences of taxa from datafile                        
// Then transform DNA sequences into numerical sequences (ppBase)   

// Compare two sequence names, treating a space or the end of the string
// as the end of the name. Needed because getSequenceData pads psname
// entries with spaces, so strcmp cannot be used on them directly.
int same_seq_name(const char *a, const char *b) {
  while (*a != 0 && *a != ' ' && *b != 0 && *b != ' ') {
    if (*a != *b) return 0;
    a++; b++;
  }
  return (*a == 0 || *a == ' ') && (*b == 0 || *b == ' ');
}

naym* getSequenceData(FILE *infile) {
  
  int i, k;
  char **ppDNASeq;
  int charLen;
  
  charLen=(nsite+1)*sizeof(char);
  ppDNASeq=(char**)malloc(nseq*sizeof(char*));
  if(ppDNASeq==NULL){
    return NULL;
  }
  
  for (i=0; i<nseq; i++)
    {
      ppDNASeq[i]=(char*)malloc(charLen);
      if(ppDNASeq[i]==NULL){
	return NULL;
      }
    }
  
  psname = (naym*)malloc(nseq*sizeof(naym));
  
  if(psname==NULL){
    return NULL;
  }

  for (i=0; i<nseq; i++)
    {
      if(fscanf(infile, "%s %s", psname[i],ppDNASeq[i]) != 2){
        printf("There was a problem reading the sequence data from the file data.phy ... exiting.\n\n");
        exit(1);
      }
      //if (verbose==1) printf("Read sequence %d with name %s\n",i,psname[i]);

      // No need to pad or load names, but done for using with other code.
      for (k=0; k<20; k++) if (psname[i][k] == 0) psname[i][k] = ' ';
 
    }
  
  transfer(ppDNASeq,ppBase,0,nsite);
  
  fclose(infile);
  
  for (i=0; i<nseq; i++) {
    free(ppDNASeq[i]);
  }
  
  free(ppDNASeq);
  
  return psname;
  
}




/*** Search all unique sites and count them - counts  are   ***/
/*** stored in the array site_counter.                      ***/

void unique_sites(){

  int i,j,m=0;
  int num_no_gaps=0;

  /*** Code efficiency from Claude, August 2026 ***/
  /* Hash-table-based rewrite of the original O(nsite^2 * nseq) all-pairs
   * comparison (each site compared against every later site to find its
   * duplicates). That algorithm is exact and was left untouched in its
   * *decision logic* -- this rewrite reaches identical site_counter[]/
   * ppBase_unique[][]/num_unique values, in the same first-occurrence
   * order, for both include_gaps settings (verified by construction: see
   * below) -- it just avoids re-scanning every prior site column by
   * hashing on a site's first encounter instead. Complexity drops from
   * O(nsite^2 * nseq) to O(nsite * nseq)
   */

  long table_size = 1;
  while (table_size < 2L*nsite + 16) table_size <<= 1;   /* next power of 2, load factor <= 0.5 */
  long table_mask = table_size - 1;
  long *slot_group = (long*)malloc(table_size*sizeof(long));
  for (long t=0; t<table_size; t++) slot_group[t] = -1;

  for (j=0; j<nsite; j++) {

    if (include_gaps==0 && ppBase[0][j]==4) continue;  /* never seeds, never matches -- always excluded */

    /* FNV-1a style hash over this site's column of base codes */
    unsigned long long h = 1469598103934665603ULL;
    for (i=0; i<nseq; i++) {
      h ^= (unsigned long)(unsigned int)ppBase[i][j];
      h *= 1099511628211ULL;
    }

    long slot = (long)(h & (unsigned long)table_mask);
    int found = -1;
    while (slot_group[slot] != -1) {
      int cand = (int)slot_group[slot];
      int match = 1;
      for (i=0; i<nseq; i++) {
        if (ppBase_unique[i][cand] != ppBase[i][j]) { match = 0; break; }
      }
      if (match) { found = cand; break; }
      slot = (slot + 1) & table_mask;
    }

    if (found >= 0) {

      site_counter[found]++;

    } else {

      site_counter[m]=1;
      for (i=0; i<nseq; i++) ppBase_unique[i][m] = ppBase[i][j];
      slot_group[slot] = m;

      m++;

    }

  }

  free(slot_group);

  num_unique=m;
  printf("Finding unique site patterns ... ");
  printf("the number of unique site patterns is %d\n\n",num_unique);
  for (i=0; i<num_unique; i++) num_no_gaps += site_counter[i];
  //for (i=0; i<num_unique; i++) printf("%d ",site_counter[i]);
  //printf("\n\n");
  if (include_gaps == 0) printf("Sites with a gap in at least one taxa have been excluded - the total number of sites used in the computations is %d\n\n",num_no_gaps);

}

void init_pattern_index(void) {
    pattern_index[0][0]=0;   pattern_index[0][1]=1;   pattern_index[0][2]=1;   pattern_index[0][3]=1;
    pattern_index[0][4]=2;   pattern_index[0][5]=7;   pattern_index[0][6]=12;  pattern_index[0][7]=12;
    pattern_index[0][8]=2;   pattern_index[0][9]=12;  pattern_index[0][10]=7;  pattern_index[0][11]=12;
    pattern_index[0][12]=2;  pattern_index[0][13]=12; pattern_index[0][14]=12; pattern_index[0][15]=7;

    pattern_index[1][0]=3;   pattern_index[1][1]=5;   pattern_index[1][2]=8;   pattern_index[1][3]=8;
    pattern_index[1][4]=6;   pattern_index[1][5]=4;   pattern_index[1][6]=10;  pattern_index[1][7]=10;
    pattern_index[1][8]=9;   pattern_index[1][9]=11;  pattern_index[1][10]=13; pattern_index[1][11]=14;
    pattern_index[1][12]=9;  pattern_index[1][13]=11; pattern_index[1][14]=14; pattern_index[1][15]=13;

    pattern_index[2][0]=3;   pattern_index[2][1]=8;   pattern_index[2][2]=5;   pattern_index[2][3]=8;
    pattern_index[2][4]=9;   pattern_index[2][5]=13;  pattern_index[2][6]=11;  pattern_index[2][7]=14;
    pattern_index[2][8]=6;   pattern_index[2][9]=10;  pattern_index[2][10]=4;  pattern_index[2][11]=10;
    pattern_index[2][12]=9;  pattern_index[2][13]=14; pattern_index[2][14]=11; pattern_index[2][15]=13;

    pattern_index[3][0]=3;   pattern_index[3][1]=8;   pattern_index[3][2]=8;   pattern_index[3][3]=5;
    pattern_index[3][4]=9;   pattern_index[3][5]=13;  pattern_index[3][6]=14;  pattern_index[3][7]=11;
    pattern_index[3][8]=9;   pattern_index[3][9]=14;  pattern_index[3][10]=13; pattern_index[3][11]=11;
    pattern_index[3][12]=6;  pattern_index[3][13]=10; pattern_index[3][14]=10; pattern_index[3][15]=4;

    pattern_index[4][0]=4;   pattern_index[4][1]=6;   pattern_index[4][2]=10;  pattern_index[4][3]=10;
    pattern_index[4][4]=5;   pattern_index[4][5]=3;   pattern_index[4][6]=8;   pattern_index[4][7]=8;
    pattern_index[4][8]=11;  pattern_index[4][9]=9;   pattern_index[4][10]=13; pattern_index[4][11]=14;
    pattern_index[4][12]=11; pattern_index[4][13]=9;  pattern_index[4][14]=14; pattern_index[4][15]=13;

    pattern_index[5][0]=7;   pattern_index[5][1]=2;   pattern_index[5][2]=12;  pattern_index[5][3]=12;
    pattern_index[5][4]=1;   pattern_index[5][5]=0;   pattern_index[5][6]=1;   pattern_index[5][7]=1;
    pattern_index[5][8]=12;  pattern_index[5][9]=2;   pattern_index[5][10]=7;  pattern_index[5][11]=12;
    pattern_index[5][12]=12; pattern_index[5][13]=2;  pattern_index[5][14]=12; pattern_index[5][15]=7;

    pattern_index[6][0]=13;  pattern_index[6][1]=9;   pattern_index[6][2]=11;  pattern_index[6][3]=14;
    pattern_index[6][4]=8;   pattern_index[6][5]=3;   pattern_index[6][6]=5;   pattern_index[6][7]=8;
    pattern_index[6][8]=10;  pattern_index[6][9]=6;   pattern_index[6][10]=4;  pattern_index[6][11]=10;
    pattern_index[6][12]=14; pattern_index[6][13]=9;  pattern_index[6][14]=11; pattern_index[6][15]=13;

    pattern_index[7][0]=13;  pattern_index[7][1]=9;   pattern_index[7][2]=14;  pattern_index[7][3]=11;
    pattern_index[7][4]=8;   pattern_index[7][5]=3;   pattern_index[7][6]=8;   pattern_index[7][7]=5;
    pattern_index[7][8]=14;  pattern_index[7][9]=9;   pattern_index[7][10]=13; pattern_index[7][11]=11;
    pattern_index[7][12]=10; pattern_index[7][13]=6;  pattern_index[7][14]=10; pattern_index[7][15]=4;

    pattern_index[8][0]=4;   pattern_index[8][1]=10;  pattern_index[8][2]=6;   pattern_index[8][3]=10;
    pattern_index[8][4]=11;  pattern_index[8][5]=13;  pattern_index[8][6]=9;   pattern_index[8][7]=14;
    pattern_index[8][8]=5;   pattern_index[8][9]=8;   pattern_index[8][10]=3;  pattern_index[8][11]=8;
    pattern_index[8][12]=11; pattern_index[8][13]=14; pattern_index[8][14]=9;  pattern_index[8][15]=13;

    pattern_index[9][0]=13;  pattern_index[9][1]=11;  pattern_index[9][2]=9;   pattern_index[9][3]=14;
    pattern_index[9][4]=10;  pattern_index[9][5]=4;   pattern_index[9][6]=6;   pattern_index[9][7]=10;
    pattern_index[9][8]=8;   pattern_index[9][9]=5;   pattern_index[9][10]=3;  pattern_index[9][11]=8;
    pattern_index[9][12]=14; pattern_index[9][13]=11; pattern_index[9][14]=9;  pattern_index[9][15]=13;

    pattern_index[10][0]=7;  pattern_index[10][1]=12; pattern_index[10][2]=2;  pattern_index[10][3]=12;
    pattern_index[10][4]=12; pattern_index[10][5]=7;  pattern_index[10][6]=2;  pattern_index[10][7]=12;
    pattern_index[10][8]=1;  pattern_index[10][9]=1;  pattern_index[10][10]=0; pattern_index[10][11]=1;
    pattern_index[10][12]=12;pattern_index[10][13]=12;pattern_index[10][14]=2; pattern_index[10][15]=7;

    pattern_index[11][0]=13; pattern_index[11][1]=14; pattern_index[11][2]=9;  pattern_index[11][3]=11;
    pattern_index[11][4]=14; pattern_index[11][5]=13; pattern_index[11][6]=9;  pattern_index[11][7]=11;
    pattern_index[11][8]=8;  pattern_index[11][9]=8;  pattern_index[11][10]=3; pattern_index[11][11]=5;
    pattern_index[11][12]=10;pattern_index[11][13]=10;pattern_index[11][14]=6; pattern_index[11][15]=4;

    pattern_index[12][0]=4;  pattern_index[12][1]=10; pattern_index[12][2]=10; pattern_index[12][3]=6;
    pattern_index[12][4]=11; pattern_index[12][5]=13; pattern_index[12][6]=14; pattern_index[12][7]=9;
    pattern_index[12][8]=11; pattern_index[12][9]=14; pattern_index[12][10]=13;pattern_index[12][11]=9;
    pattern_index[12][12]=5; pattern_index[12][13]=8; pattern_index[12][14]=8; pattern_index[12][15]=3;

    pattern_index[13][0]=13; pattern_index[13][1]=11; pattern_index[13][2]=14; pattern_index[13][3]=9;
    pattern_index[13][4]=10; pattern_index[13][5]=4;  pattern_index[13][6]=10; pattern_index[13][7]=6;
    pattern_index[13][8]=14; pattern_index[13][9]=11; pattern_index[13][10]=13;pattern_index[13][11]=9;
    pattern_index[13][12]=8; pattern_index[13][13]=5; pattern_index[13][14]=8; pattern_index[13][15]=3;

    pattern_index[14][0]=13; pattern_index[14][1]=14; pattern_index[14][2]=11; pattern_index[14][3]=9;
    pattern_index[14][4]=14; pattern_index[14][5]=13; pattern_index[14][6]=11; pattern_index[14][7]=9;
    pattern_index[14][8]=10; pattern_index[14][9]=10; pattern_index[14][10]=4; pattern_index[14][11]=6;
    pattern_index[14][12]=8; pattern_index[14][13]=8; pattern_index[14][14]=5; pattern_index[14][15]=3;

    pattern_index[15][0]=7;  pattern_index[15][1]=12; pattern_index[15][2]=12; pattern_index[15][3]=2;
    pattern_index[15][4]=12; pattern_index[15][5]=7;  pattern_index[15][6]=12; pattern_index[15][7]=2;
    pattern_index[15][8]=12; pattern_index[15][9]=12; pattern_index[15][10]=7; pattern_index[15][11]=2;
    pattern_index[15][12]=1; pattern_index[15][13]=1; pattern_index[15][14]=1; pattern_index[15][15]=0;
}


int main(int argc, char *argv[]) {

  FILE *set, *tf, *data, *res;
  char keyword;
  int memcheck, i, j, k, l;
  int i2, j2, k2, l2;
  int boottype, nboot, ntquarts;
  int table_initialized = 0;
  naym tempname1, tempname2;
  naym *map_seqname = NULL;   // sequence names from the settings file assignments
  int *map_species = NULL;    // species index for each entry of map_seqname
  int *map_used = NULL;       // flags settings entries already matched to data.phy
  float temp_rate;
  double totaltime, beta_init;

  /* --- defaults --- */
  const char *settings_file   = "settings";
  const char *data_file       = "data.phy";
  const char *tree_file       = "treefile.tre";
  const char *out_file        = "outtree.tre";
  const char *picltrees_file  = "picltrees.tre";
  const char *results_file    = "results";
  const char *bootdata_file   = "boots.dat";


  /* --- override via argv --- */
  if (argc > 1) settings_file  = argv[1];
  if (argc > 2) data_file      = argv[2];
  if (argc >3) tree_file       = argv[3];
  if (argc > 4) out_file       = argv[4];
  if (argc > 5) picltrees_file = argv[5];
  if (argc > 6) results_file   = argv[6];
   if (argc > 7) bootdata_file = argv[7];

  if(argc >1) {
    setvbuf(stdout, NULL, _IOLBF, 0);   /* line-buffered, even on a pipe */
    setvbuf(stderr, NULL, _IONBF, 0);   /* unbuffered for errors        */

    printf("Using files:\n");
    printf("  settings:    %s\n", settings_file);
    printf("  data:        %s\n", data_file);
    printf("  treefile:    %s\n", tree_file);
    printf("  outtree:     %s\n", out_file);
    printf("  picltrees:   %s\n", picltrees_file);
    printf("  results:     %s\n\n", results_file);
    printf("  bootdata:    %s\n\n", bootdata_file);
 }

  // print opening message
  printf("\n\n----------------------------------------------------------------------------\n");
  printf("\t PICL: Phylogenetic Inference with Composite Likelihood \n");
  printf("\t\t\t     Version 1.2  \n");
  printf("\t\t\t   October 17, 2026 \n");
  printf("----------------------------------------------------------------------------\n\n");

  /*******************************************/
  /* Read in settings from the settings file */
  /*******************************************/

  // read from settings file
  set = fopen(settings_file,"r");
  if (set==NULL) { printf("Settings file not found ... exiting.\n\n"); exit(1);}
  char mykeyword[15];
  
  printf("Model parameters are now being read from the settings file. Please check the listing\n");
  printf("below to make sure that your settings are being read in as intended:\n\n");

  /* Model */
  if (fscanf(set,"%s %d",mykeyword,&model) != 2) {
	printf("\nFailed to read from settings file -- error on line 1\n");
        printf("Exiting.\n\n");
    	exit(1);
  }
  if (strcmp(mykeyword,"Model:") == 0)	printf("%s %d\n",mykeyword,model);
  else {
	printf("\nThe settings file should have contained the keyword 'Model' followed by a colon on line 1.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
	printf("Exiting.\n\n");
	exit(1);
  }

  /* Gaps */
  if (fscanf(set,"%s %d",mykeyword,&include_gaps) != 2) { 
        printf("\nFailed to read from settings file -- error on line 2\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Gaps:") == 0)  printf("%s %d\n",mykeyword,include_gaps);
  else {
        printf("\nThe settings file should have contained the keyword 'Gaps' followed by a colon on line 2.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }


  /* Boot_type */
  if (fscanf(set,"%s %d",mykeyword,&boottype) != 2) {
        printf("\nFailed to read from settings file -- error on line 3\n");
	printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Boot_type:") == 0)  printf("%s %d\n",mykeyword,boottype);
  else {
        printf("\nThe settings file should have contained the keyword 'Boot_type' followed by a colon on line 3.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");    
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Nboot */
  if (fscanf(set,"%s %d",mykeyword,&nboot) != 2) {
        printf("\nFailed to read from settings file -- error on line 4\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Nboot:") == 0)  printf("%s %d\n",mykeyword,nboot);
  else {
        printf("\nThe settings file should have contained the keyword 'Nboot' followed by a colon on line 4.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");    
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Theta */
  if (fscanf(set,"%s %f",mykeyword,&theta) != 2) {
        printf("\nFailed to read from settings file -- error on line 5\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Theta:") == 0)  printf("%s %f\n",mykeyword,theta);
  else {
        printf("\nThe settings file should have contained the keyword 'Theta' followed by a colon on line 5.\n");            
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  /**** Now set theta = 4Nmu to theta/2, which is what this code uses throughout ****/
  theta = theta/2.0;

  /* Alpha */
  if (fscanf(set,"%s %f",mykeyword,&alpha) != 2) {
        printf("\nFailed to read from settings file -- error on line 6\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Alpha:") == 0)  printf("%s %f\n",mykeyword,alpha);
  else {
        printf("\nThe settings file should have contained the keyword 'Alpha' followed by a colon on line 6.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Lambda */
  if (fscanf(set,"%s %f",mykeyword,&lambda) != 2) {
        printf("\nFailed to read from settings file -- error on line 7\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Lambda:") == 0)  printf("%s %f\n",mykeyword,lambda);
  else {
        printf("\nThe settings file should have contained the keyword 'Lambda' followed by a colon on line 7.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
	
  /* Rate_param */
  if (fscanf(set,"%s %f",mykeyword,&ratepar) != 2) {
        printf("\nFailed to read from settings file -- error on line 8\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Rate_param:") == 0)  printf("%s %f\n",mykeyword,ratepar);
  else {
        printf("\nThe settings file should have contained the keyword 'Rate_param' followed by a colon on line 8.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Random_tree */
  if (fscanf(set,"%s %d",mykeyword,&random_tree) != 2) {
        printf("\nFailed to read from settings file -- error on line 9\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Random_tree:") == 0)  printf("%s %d\n",mykeyword,random_tree);
  else {
        printf("\nThe settings file should have contained the keyword 'Random_tree' followed by a colon on line 9.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");                                                                                                          
        exit(1);
  }

  /* Opt_bl */
  if (fscanf(set,"%s %d",mykeyword,&anneal_bl) != 2) {
        printf("\nFailed to read from settings file -- error on line 10\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Opt_bl:") == 0)  printf("%s %d\n",mykeyword,anneal_bl);
  else {
        printf("\nThe settings file should have contained the keyword 'Opt_bl' followed by a colon on line 10.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* User_bl */
  if (fscanf(set,"%s %d",mykeyword,&user_bl) != 2) {
        printf("\nFailed to read from settings file -- error on line 11\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"User_bl:") == 0)  printf("%s %d\n",mykeyword,user_bl);
  else {
        printf("\nThe settings file should have contained the keyword 'User_bl' followed by a colon on line 11.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Num_opt */
  if (fscanf(set,"%s %d",mykeyword,&max_it_bl) != 2) {
        printf("\nFailed to read from settings file -- error on line 12\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Num_opt:") == 0)  printf("%s %d\n",mykeyword,max_it_bl);
  else {
        printf("\nThe settings file should have contained the keyword 'Num_opt' followed by a colon on line 12.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Seed1 */
  if (fscanf(set,"%s %d",mykeyword,&seedj) != 2) {  
        printf("\nFailed to read from settings file -- error on line 13\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Seed1:") == 0)  printf("%s %d\n",mykeyword,seedj);      
  else {
        printf("\nThe settings file should have contained the keyword 'Seed1' followed by a colon on line 13.\n");   
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Seed2 */
  if (fscanf(set,"%s %d",mykeyword,&seedk) != 2) {      
        printf("\nFailed to read from settings file -- error on line 14\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Seed2:") == 0)  printf("%s %d\n",mykeyword,seedk);            
  else {
        printf("\nThe settings file should have contained the keyword 'Seed2' followed by a colon on line 14.\n");     
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Num_cat */
  if (fscanf(set,"%s %d",mykeyword,&ncat) != 2) {      
        printf("\nFailed to read from settings file -- error on line 15\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Num_cat:") == 0)  printf("%s %d\n",mykeyword,ncat);            
  else {
        printf("\nThe settings file should have contained the keyword 'Num_cat' followed by a colon on line 15.\n");     
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Tree_search */
  if (fscanf(set,"%s %d",mykeyword,&anneal) != 2) {      
        printf("\nFailed to read from settings file -- error on line 16\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Tree_search:") == 0)  printf("%s %d\n",mykeyword,anneal);            
  else {
        printf("\nThe settings file should have contained the keyword 'Tree_search' followed by a colon on line 16.\n");     
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  /* Max_iter */
  if (fscanf(set,"%s %d",mykeyword,&max_it) != 2) {
        printf("\nFailed to read from settings file -- error on line 17\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Max_iter:") == 0)  printf("%s %d\n",mykeyword,max_it);
  else {
        printf("\nThe settings file should have contained the keyword 'Max_iter' followed by a colon on line 17.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
 
  /* Multi_iter */                                                                    
  if (fscanf(set,"%s %d",mykeyword,&mult_iter) != 2) {
        printf("\nFailed to read from settings file -- error on line 18\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Multi_iter:") == 0)  printf("%s %d\n",mykeyword,mult_iter);
  else {
        printf("\nThe settings file should have contained the keyword 'Multi_iter' followed by a colon on line 18.\n");
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  //mult_iter = 3;

  /* prob_bound */
  if (fscanf(set,"%s %lf",mykeyword,&prob_bound) != 2) {
        printf("\nFailed to read from settings file -- error on line 19\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Prob_bound:") == 0)  printf("%s %f\n",mykeyword,prob_bound);
  else {
        printf("\nThe settings file should have contained the keyword 'Prob_bound' followed by a colon on line 19.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  //prob_bound = 0.05;

  /* Test_incr */
  if (fscanf(set,"%s %d",mykeyword,&test_increment) != 2) {
        printf("\nFailed to read from settings file -- error on line 20\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Test_incr:") == 0)  printf("%s %d\n",mykeyword,test_increment);
  else {
        printf("\nThe settings file should have contained the keyword 'Test_incr' followed by a colon on line 20.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  //test_increment = 250; 

  /* Opt_slope */
  if (fscanf(set,"%s %lf",mykeyword,&b1opt) != 2) {
        printf("\nFailed to read from settings file -- error on line 21\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Opt_slope:") == 0)  printf("%s %f\n",mykeyword,b1opt);
  else {
        printf("\nThe settings file should have contained the keyword 'Opt_slope' followed by a colon on line 21.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  //b1opt = -0.01;
  
  /* Beta */
  if (fscanf(set,"%s %f",mykeyword,&beta) != 2) {
        printf("\nFailed to read from settings file -- error on line 22\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Beta:") == 0)  printf("%s %f\n",mykeyword,beta);
  else {
        printf("\nThe settings file should have contained the keyword 'Beta' followed by a colon on line 22.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  beta_init = beta;  /* used to initialize bootstrapping search */

  /* Verbose */
  if (fscanf(set,"%s %d",mykeyword,&verbose) != 2) {
        printf("\nFailed to read from settings file -- error on line 23\n");
        printf("Exiting.\n\n");
        exit(1);
  }
  if (strcmp(mykeyword,"Verbose:") == 0)  printf("%s %d\n",mykeyword,verbose);
  else {
        printf("\nThe settings file should have contained the keyword 'Verbose' followed by a colon on line 23.\n"); 
        printf("Please correct your settings file. Note that the keyword order must match that in the example file.\n\n");
        printf("Exiting.\n\n");
        exit(1);
  }

  if (model == 5)theta = NAN;  	  
  if (b1opt>0) printf("The target annealing slope (b1opt) must be negative; adaptive annealing is disabled.\n\n");

  /*****************************/
  /* End read model settings   */
  /*****************************/


  if(fscanf(set,"%d",&ntaxa) != 1 || ntaxa <= 0){
    printf("Invalid ntaxa");
    exit(1);
  }

  setall((long)seedj,(long)seedk); /* Set seeds for random number generator */
  nquarts = binomial(ntaxa,4);
  printf("\nThere are %d total species. This corresponds to %d species-level quartets.\n",ntaxa,nquarts);
  data = fopen(data_file,"r");
  if( fscanf(data,"%d %d\n",&nseq,&nsite) != 2 || nseq <= 0 || nsite <= 0){
    printf("Invalid number of sequences or sites");
    exit(1);
  }
  ntquarts = binomial(nseq,4);
  printf("There are %d total sequences. This corresponds to %d sequence-level quartets.\n\n",nseq,ntquarts);

  // allocate memory
  memcheck = MemAlloc();
  for (i=0; i<2*ntaxa; i++) {
    filled_ind[i] = 0;
  }
  for (i=0;i<(ntaxa+1);i++){
      ppTwoRow_temp[0][i] = 0;
      ppTwoRow_temp[1][i] = 0;
      seq_counter[i] = 0;
  }

  if(ntaxa <= 0 || ntaxa >100000) {
    printf("Too many taxa or negative taxa");
    exit(1);
  }
  else {
    taxname = (naym*)malloc(ntaxa*sizeof(naym));
  }
  
  for (i=0; i<ntaxa; i++) {
    if(fscanf(set,"%s",taxname[i]) != 1){
      printf("The number of species names specified in the settings file is incorrect ... exiting.\n\n");
      exit(1);
    }
  } 


  // Read the species/sequence-name pairs from the settings file.
  // These are stored by name only; the mapping to rows of data.phy is
  // done after the alignment is read, so data.phy may be in any order.
  map_seqname = (naym*)malloc(nseq*sizeof(naym));
  map_species = (int*)malloc(nseq*sizeof(int));
  if (map_seqname==NULL || map_species==NULL) {
    printf("Can't memalloc species/sequence map ... exiting.\n\n");
    exit(1);
  }

  for (i=0; i<nseq; i++) {
        if( fscanf(set,"%s %s",tempname1,tempname2) != 2){
          printf("The number of taxon assignments in the settings file is incorrect (expected %d, one per sequence in %s) ... exiting.\n\n",nseq,data_file);
          exit(1);
        }

        // find the species this sequence belongs to
        for (j=0; j<ntaxa; j++) if (strcmp(tempname1,taxname[j]) == 0) break;
        if (j == ntaxa) {
          printf("Species name '%s' in the taxon assignments of the settings file is not among the %d species names listed ... exiting.\n\n",tempname1,ntaxa);
          exit(1);
        }

        // make sure the sequence name is not assigned twice
        for (k=0; k<i; k++) {
          if (strcmp(tempname2,map_seqname[k]) == 0) {
            printf("Sequence name '%s' is assigned more than once in the settings file ... exiting.\n\n",tempname2);
            exit(1);
          }
        }

        strcpy(map_seqname[i],tempname2);
        map_species[i] = j;
  }

  if (model==1) printf("Estimation will be carried out using the standard CIS/multilocus model (MSC-JC69)\n");
  else if (model==2) printf("Estimation will be carried out using the CIS/multlilocus rate variation model (MSC-JC69 + G)\n");
  else if (model==3) printf("Estimation will be carried out using the SNP model\n");
  else if (model==4) printf("Estimation will be carried out using the JC69 model for gene trees\n");
  else if (model==5) printf("Estimation will be carried out using the CIS/multilocus model (MSC-JC69) with variable population size\n");
  else { printf("Model must be 1, 2, 3, 4 or 5. Exiting.\n\n"); exit(1); }
 
  /******* TO DO: need to free extra memory from storing taxon map *******/

  // set up look up table to deal with ambiguity codes efficiently

  if (!table_initialized) {
    for (int code = 0; code < 15; code++) {
        BaseWeights bw = get_base_weights(code);
        for (int k = 0; k < 4; k++)
            base_weight_table[code][k] = bw.weight[k];
    }
    table_initialized = 1;
  }

  // read in data
  getSequenceData(data);
  // fclose(data);
  printf("\nSequence alignment with %d lineages and %d sites has been read successfully\n",nseq,nsite);

  // Map each sequence (row of data.phy) to its species by matching names
  // against the assignments read from the settings file.
  map_used = (int*)calloc(nseq,sizeof(int));
  if (map_used==NULL) {
    printf("Can't memalloc map_used ... exiting.\n\n");
    exit(1);
  }

  for (i=0; i<nseq; i++) {
        for (k=0; k<nseq; k++) if (same_seq_name(psname[i],map_seqname[k])) break;
        if (k == nseq) {
          // print the name only up to the padding added in getSequenceData
          for (l=0; psname[i][l] != ' ' && psname[i][l] != 0; l++);
          printf("Sequence '%.*s' in %s has no species assignment in the settings file ... exiting.\n\n",l,psname[i],data_file);
          exit(1);
        }
        if (map_used[k]) {
          printf("Sequence name '%s' appears more than once in %s ... exiting.\n\n",map_seqname[k],data_file);
          exit(1);
        }
        map_used[k] = 1;

        j = map_species[k];
        if (seq_counter[j] >= 100) {
          printf("Species %s has more than 100 sequences, which exceeds the space allocated in ppSp_assign ... exiting.\n\n",taxname[j]);
          exit(1);
        }
        ppSp_assign[j][seq_counter[j]] = i;
        seq_counter[j]++;
  }

  // every species must be represented by at least one sequence
  for (j=0; j<ntaxa; j++) {
        if (seq_counter[j] == 0) {
          printf("Species %s has no sequences assigned to it ... exiting.\n\n",taxname[j]);
          exit(1);
        }
  }

  // store the sequence names belonging to each species
  assignment_vector = (naym**)malloc(ntaxa*sizeof(naym*));
  if (assignment_vector==NULL) {
    printf("Can't memalloc assignment_vector ... exiting.\n\n");
    exit(1);
  }
  for (j=0; j<ntaxa; j++) {
        assignment_vector[j] = (naym*)malloc(seq_counter[j]*sizeof(naym));
        if (assignment_vector[j]==NULL) {
          printf("Can't memalloc assignment_vector[%d] ... exiting.\n\n",j);
          exit(1);
        }
  }
  for (k=0; k<nseq; k++) {
        // find the row in data.phy that this settings entry was matched to
        j = map_species[k];
        for (i=0; i<seq_counter[j]; i++) {
          if (same_seq_name(psname[ppSp_assign[j][i]],map_seqname[k])) {
            strcpy(assignment_vector[j][i],map_seqname[k]);
            break;
          }
        }
  }

  free(map_seqname);
  free(map_species);
  free(map_used);

  // checks that assignments were read in correctly
  for (j=0; j<ntaxa; j++) printf("Species %d with name %s has %d sequences\n",j,taxname[j],seq_counter[j]);
  printf("\n");
  fflush(0);

  printf("The species assignments are:\n");
  for (j=0; j<ntaxa; j++) {
	printf("Species %d: %s\n",j,taxname[j]);
	for (i=0; i<seq_counter[j]; i++) {
		printf("%d (%s) ",ppSp_assign[j][i],assignment_vector[j][i]);
	}
	printf("\n");
  }
  printf("\n");
  unique_sites();

  // build lookup table for quartet site pattern counts
  init_pattern_index();

  // set a tree
  if (random_tree==0) { // read in user species tree
  	for (i=0; i<2*ntaxa; i++) TimeVec[i] = 0.0;
  	tf = fopen(tree_file,"r");
	if (tf==NULL) { printf("File treefile.tre not found. Exiting.\n\n"); exit(1);}
  	temp_rate = ReadTree(tf);
  
  	totaltime = FindTotalTime();
  	
  	CalcTimeVec(totaltime,1.0);
	
  	// Transfer from temp and re-scale with theta if needed

  	for (i=0; i<ntaxa-1; i++) {
		ppTwoRow[0][i] = ppTwoRow_temp[0][i];
	 	ppTwoRow[1][i] = ppTwoRow_temp[1][i];
  	}
	if (user_bl==0) { //do not use times from user tree; generate equally-spaced times
		compute_times(ntaxa+1);
		for (i=0; i<2*ntaxa+1; i++) TimeVec_temp[i] = TimeVec[i];
	}
	if(model != 5){  /* convert to mutation units */
	    for (i=0; i<2*ntaxa+1; i++) {
		TimeVec[i] = theta*TimeVec_temp[i];
		TimeVec_init[i] = theta*TimeVec_temp[i];
		TimeVec_temp[i] = TimeVec[i];
  	    }
	}
  	else if(model ==5) {
		for (i=0; i<2*ntaxa+1; i++) {
		TimeVec_init[i] = TimeVec[i];
		TimeVec_temp[i] = TimeVec[i];
		}
	}

  	if (verbose == 1) {
  	    if(model !=5){
  	        printf("The tree read from the input file is (coalescent units, mutation units):\n");
      	    for (i=0; i<ntaxa-1; i++) {
      		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
      		printf("%f %f\n ",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
      		}
      		printf("\n");
  	    }
  	    else if(model ==5){
  	        printf("The tree read from the input file is (mutation units):\n");
      	    for (i=0; i<ntaxa-1; i++) {
      		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
      		printf("%f \n ",TimeVec[i+ntaxa+1]);
      		}
  	    }
    		
    		
		
    }
  	printf("Tree successfully read from file\n\n");        
   }
   else if (random_tree==1) {  // generate a starting tree under the Yule model

	build_random_tree(seedj,seedk);
	/* populate temp structures */
        for (i=0; i<ntaxa-1; i++) {
                ppTwoRow_temp[0][i] = ppTwoRow[0][i];
                ppTwoRow_temp[1][i] = ppTwoRow[1][i];
        }
	for (i=0; i<ntaxa+1; i++) TimeVec[i] = 0.0;
	for (i=0; i<2*ntaxa+1; i++) {
		TimeVec_temp[i] = TimeVec[i]; 
                TimeVec_init[i] = TimeVec[i];
        }
 	if (verbose == 1 ) {
				if (model == 5){
					printf("The random tree generated is (mutation units):\n");
                	for (i=0; i<ntaxa-1; i++) {
               		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
                	printf("%f\n ",TimeVec[i+ntaxa+1]);
				}
				}
				else{
                printf("The random tree generated is (coalescent units, mutation units):\n");
                for (i=0; i<ntaxa-1; i++) {
                printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
                printf("%f %f\n ",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
				}
        	}
	printf("\n");
	}
   }
   else {  // print an error and exit
	printf("Random_tree must be either 0 or 1. Exiting.\n\n");
	exit(1);
   }
   fflush(0);

  // make the IND matrix that will be used for composite likelihood computation
  //make_indmat();
  for (i=1; i<2*ntaxa; i++) {
	parents[i] = find_parent_scan(i);
 	parents_temp[i] = parents[i];
  }
  make_indmat();

  /* test -- getting estimates of branch lengths via computation -- feature not yet implemented */
  //AllMyChildrenLeft(8);
  //AllMyChildrenRight(8);
  //AllMyChildren(8);
  //GetMomentEstimate(8);
  //for (i=0; i<3; i++) for (j=0;j<ntaxa;j++) printf("%d ",ppNodeChildren[i][j]);
  /* end test */

  /* Compute the composite likelihood for the current tree */
  if (model == 5) {
	  ComputeAandS_popvar(alpha, lambda);
  }
  else {
  	ComputeAandS(theta);
  }
  for (i=0; i<nquarts+1; i++) qvec[i]=0; 
  if (model == 1) curr_anneal_lik = GetCompLik();
  else if (model == 2) curr_anneal_lik = GetCompLik_ratevar();
  else if (model == 3) curr_anneal_lik = GetCompLik_msnp();
  else if (model == 4) curr_anneal_lik = GetCompLik_genetree(); 
  else if (model == 5) curr_anneal_lik = GetCompLik_popvar(); 
  printf("The composite likelihood of the initial tree is %f\n\n",curr_anneal_lik);


  /*************  Set-up complete ***************/


  /*************  Run the appropriate analysis **************/

  // Evalute current tree only
  // if anneal = 0
  if (anneal==0) {
	if (anneal_bl==0) { printf("Analysis complete -- exiting.\n\n"); exit(1); }
	else if (anneal_bl==1) { // uphill bl search
		printf("Optimizing branch lengths using uphill search ...");
		if (model == 1) bl_uphill_full();
        	else if (model == 2) {
                	bl_uphill_ratevar();                                                                                               
                	OptAlpha();
                	bl_uphill_ratevar();
        	}
        	else if (model == 3) bl_uphill_msnp();
		else if (model == 4) bl_uphill_genetree();
		else if (model == 5) bl_uphill_full_popvar();
	}
	else if (anneal_bl==2) { // simulated annealing bl search
		printf("Optimizing branch lengths using simulating annealing ...");
		if (model == 1) bl_anneal_full();
        	else if (model == 2) {
                	bl_anneal_ratevar();
                	OptAlpha();
                	bl_anneal_ratevar();
        	}
        	else if (model == 3) bl_anneal_msnp();
		else if (model == 4) bl_anneal_genetree();
		else if (model == 5) bl_anneal_full_popvar();
	}
	else if (anneal_bl==3) { // numerical optimization of branch lengths
		printf("Numerical optimization of branch lengths is not yet implemented. Exiting.\n\n"); exit(1);
	}
	else { printf("Opt_bl must be 0, 1, 2, or 3. Exiting.\n\n"); exit(1); }

        if (verbose == 1 ) {
			if (model == 5){
				printf("The tree after branch length optimization is (mutation units):\n");
    			for (i=0; i<ntaxa-1; i++) {
        				printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        				printf("%f \n",TimeVec[i+ntaxa+1]);
					}
			}
			else{
	    		printf("The tree after branch length optimization is (coalescent units, mutation units):\n");
	    		for (i=0; i<ntaxa-1; i++) {
	        		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
	        		printf("%f %f\n",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
				}
			}
		}

	if (model == 1) printf(" the composite likelihood of the tree is %f\n",GetCompLik());
  	else if (model == 2) printf(" the composite likelihood of the tree is %f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);
  	else if (model == 3) printf(" the composite likelihood of the tree is %f\n",GetCompLik_msnp());
	else if (model == 4) printf(" the composite likelihood of the tree is %f\n",GetCompLik_genetree());
	else if (model == 5) printf(" the composite likelihood of the tree is %f\n",GetCompLik_popvar());
  	printf("\n\n");
  	
	// re-order ppTwoRow and write tree to outtree.tre
        int temptax;
        for (i=0; i<ntaxa; i++) {
                if (ppTwoRow[0][i]>ppTwoRow[1][i]) {
                temptax = ppTwoRow[0][i];
                ppTwoRow[0][i] = ppTwoRow[1][i];
                ppTwoRow[1][i] = temptax;
                }
        }

        out = fopen(out_file,"w");
        fprintf(out,"Current tree in mutation units: \n");
        write_species_tree_out(ntaxa+1,ntaxa+1);
	fprintf(out,";\n\n");
	if (model != 5) {
		fprintf(out,"Current tree in coalescent units: \n");
        	for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta;
        	write_species_tree_out(ntaxa+1,ntaxa+1);
        	fprintf(out,";\n\n");
        	for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]*theta;
	}
	printf("Tree has been written to file outtree.tre\n\n");
  }

  // Uphill tree search 
  // if anneal = 1
  else if (anneal==1) {

	printf("Uphill tree searching is not yet implemented -- exiting.\n");

  }

  // Simulated annealing tree search
  // if anneal = 2
  else if (anneal==2) {
	ci = ntaxa;
	printf("Tree is being estimated with simulated annealing ...");
        fflush(0);

	// copy starting tree to besttree before beginning annealing
	for (i=0; i<ntaxa; i++) {
                ppTwoRow_best[0][i] = ppTwoRow[0][i];
                ppTwoRow_best[1][i] = ppTwoRow[1][i];
        }
        for (i=0; i<2*ntaxa+1; i++) TimeVec_best[i] = TimeVec[i];

	// now start annealing
	if (model == 1) {
		//bl_uphill_full();
		anneal_full();
		if (anneal_bl==1) bl_uphill_full();
		else if (anneal_bl==2) bl_anneal_full();
	} 
	else if (model == 2) {
		bl_uphill_ratevar();
		anneal_ratevar();
		if (anneal_bl==1) bl_uphill_ratevar();
		else if (anneal_bl==2) bl_anneal_ratevar();
	}
	else if (model == 3) {
		bl_uphill_msnp();
		anneal_msnp();
		if (anneal_bl==1) bl_uphill_msnp();
		else if (anneal_bl==2) bl_anneal_msnp();
	}
	else if (model == 4) {
		bl_uphill_genetree();
		anneal_genetree();
		if (anneal_bl==1) bl_uphill_genetree();
		else if (anneal_bl==2) bl_anneal_genetree();
	}
	else if (model == 5) {
		//bl_uphill_full();
		anneal_full_popvar();
		if (anneal_bl==1) bl_uphill_full_popvar();
		else if (anneal_bl==2) bl_anneal_full_popvar();
	}

	//printf(" done\n\n");
	if (verbose==1) {
		printf("After search, tree is \n");
  		for (i=0; i<ntaxa-1; i++) {
        		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
				if (model == 5) {
					printf("%f \n",TimeVec[i+ntaxa+1]);
				}
				else {
					printf("%f %f\n",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
				}
        		
  		}
		if (model == 1) printf("The composite likelihood of this tree is %f\n",GetCompLik());
                else if (model == 2) printf("The composite likelihood of this tree is %f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);
                else if (model == 3) printf("The composite likelihood of this tree is %f\n",GetCompLik_msnp());
				else if (model == 5) printf("The composite likelihood of this tree is %f\n",GetCompLik_popvar());
                printf("\n\n");
	}

	// re-order ppTwoRow
  	int temptax;
  	for (i=0; i<ntaxa; i++) {   
        	if (ppTwoRow[0][i]>ppTwoRow[1][i]) {  
                temptax = ppTwoRow[0][i];
                ppTwoRow[0][i] = ppTwoRow[1][i];
                ppTwoRow[1][i] = temptax;
        	}
  	}

       	out = fopen(out_file,"w");
  	fprintf(out,"Current tree in mutation units: \n");
  	write_species_tree_out(ntaxa+1,ntaxa+1);
	fprintf(out,";\n\n");
	if (model != 5){
		fprintf(out,"Current tree in coalescent units: \n");
		for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta;
		write_species_tree_out(ntaxa+1,ntaxa+1);
  		fprintf(out,";\n\n");   
		for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]*theta; 
	}
  	// copy best tree to current tree and write to file outtree.tre
  	for (i=0; i<ntaxa; i++) {
         	ppTwoRow[0][i] = ppTwoRow_best[0][i];
         	ppTwoRow[1][i] = ppTwoRow_best[1][i];
		ppTwoRow_best_orig[0][i] = ppTwoRow_best[0][i];
		ppTwoRow_best_orig[1][i] = ppTwoRow_best[1][i];
  	}
  	for (i=0; i<2*ntaxa+1; i++) TimeVec[i] = TimeVec_best[i];
	for (i=0; i<2*ntaxa+1; i++) TimeVec_best_orig[i] = TimeVec_best[i];
	// ppTwoRow was just replaced wholesale with the best tree found during
	// the search, so parents[] (kept in sync incrementally by trbldg.c
	// during the search) is now stale relative to it -- resync from
	// scratch before the GetCompLik() call below reads it via
	// find_parent(). Missing this resync silently corrupted the reported
	// "composite likelihood of the best tree" (find_parent() would return
	// parent indices for the pre-transfer tree).
	for (i=1; i<2*ntaxa; i++) {
		parents[i] = find_parent_scan(i);
		parents_temp[i] = parents[i];
	}

	// re-order ppTwoRow to write to outtree.tre

        for (i=0; i<ntaxa; i++) {
                if (ppTwoRow[0][i]>ppTwoRow[1][i]) {
                temptax = ppTwoRow[0][i];
                ppTwoRow[0][i] = ppTwoRow[1][i];
                ppTwoRow[1][i] = temptax;
                }
        }
	if (verbose == 1){
		printf("After transfer from bessttree, tree is \n");
					if (model == 5){
						for (i=0; i<ntaxa-1; i++) {
                        	printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);  
                        	printf("%f \n",TimeVec[i+ntaxa+1]);
					}
					}
					else {
						
                	for (i=0; i<ntaxa-1; i++) {
                        	printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);  
                        	printf("%f %f\n",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
					}
        	}
	}
	printf("The composite likelihood of the best tree found by the algorithm is ");
	if (model == 1) printf("%f\n",GetCompLik());
    	else if (model == 2) printf("%f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);    
    	else if (model == 3) printf("%f\n",GetCompLik_msnp());
	else if (model == 4) printf("%f\n",GetCompLik_genetree());
	else if (model == 5) printf("%f\n",GetCompLik_popvar());
        printf("\n\n");
         
   	fprintf(out,"Best tree found by the algorithm in mutation units: \n");
  	write_species_tree_out(ntaxa+1,ntaxa+1);
  	fprintf(out,";\n\n");
	if (model != 5){
		
	fprintf(out,"Best tree found by the algorithm in coalescent units: \n");
        for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta; 
        write_species_tree_out(ntaxa+1,ntaxa+1);     
	fprintf(out,";\n\n");
	}
	fclose(out);
	printf("Results have been written to file outtree.tre.\n\n");

  }
  else { printf("Tree_search must be 0, 1, or 2. Exiting.\n\n"); exit(1); }
  /* main analyses done */


  /* print branch lengths to file */
  if (verbose==1) {
  	res = fopen(results_file,"w");
	if (model == 5) {
		for (i=1; i<ntaxa; i++) fprintf(res,"%f ",TimeVec[ntaxa+i]);
  	fclose(res);
	}
	  else {
		  for (i=1; i<ntaxa; i++) fprintf(res,"%f ",TimeVec[ntaxa+i]/theta);
  			fprintf(res,"%f \n",ratepar);
  			fclose(res);
	  }
  	
  }
  /* done print branch lengths to file */

  /* print only Newick string to file */
	pt = fopen(picltrees_file,"w");
	write_species_tree_out_file(ntaxa+1,ntaxa+1);
	fprintf(pt,";\n");
	fclose(pt);
  /* done print only Newick string to file */


  /* bootstrapping */
  if (nboot>0 && boottype==1) {
	if (anneal == 2) {
		printf("Beginning bootstrap for trees\n\n");
		boot = fopen(bootdata_file,"w");
		bootstrap(nboot,beta_init);
		fclose(boot);
	}
	else printf("Tree bootstrapping can only be done following tree estimation. Please re-run with tree search set to 2.\n\n");	
  }
  else if(nboot>0 && boottype==2) {
	if (anneal_bl == 0) { printf("Please set Opt_bl = 1 to carry out bootstrapping for the branch lengths along a fixed tree -- exiting.\n\n"); exit(1); }
	else if (anneal != 0) { printf("Branch length bootstrapping can only be done on a fixed input tree. Please set Tree_search to 0 and be sure that the tree of interest is in the file treefile.tre.\n\n"); exit(1); }
	else {
		printf("Beginnning bootstrapping for branch lengths along a fixed tree\n\n");
		boot_times(nboot,bootdata_file);
	}
  }

  printf("\n\nAnalysis complete - exiting.\n\n");

  /* useful functions -- not accessible to users through settings */
    //remove_CONSTANT();
    //remove_CONSTANT_MSNP();
    //print_PHYLIP();

}
