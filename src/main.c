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
#include "boot.c"
#include "opt.c"
#include "trbldg.c"
#include "trbldg_ratevar.c"
#include "trbldg_msnp.c"
#include "trbldg_genetree.c"
#include "complik_popvar.c"
#include "treebldg_popvar.c"
#include "anneal_popvar.c"


// ntaxa = number of species
// nseq = number of individuals
int ntaxa, nspecies, nseq, nsite, nquarts, num_unique_quarts, num_unique, include_gaps, num_no_gaps, verbose, model, ncat, random_tree;
int anneal, anneal_bl, user_bl, max_it, mult_iter, num_reject, max_it_bl, test_increment, seedj, seedk;
int *parents, *parents_temp, *ppTwoRow[2], *ppTwoRow_temp[2], *ppTwoRow_best[2], *ppTwoRowQuart[2], *filled_ind, *seq_counter, *qvec, *site_counter;
int **ppBase_full, **ppBase, **ppBase_unique, **ppSp_assign, **ppNodeChildren, **ppNodeChildrenLeftQuart, **ppNodeChildrenRightQuart;
int pattern_index[16][16];
double ci, max_cl, curr_anneal_lik, b1opt, prob_bound;
float lambda, theta, beta, mu, ratepar, invpar;
double *TimeVec, *TimeVec_temp, *TimeVec_init, *TimeVec_best, *TimeVecQuart, *rvals, **ppLengthMat, **ppMatrix;
double smat[10][10],amat[12][12];
double base_weight_table[15][4]; 
FILE *out, *pt;

Link Head, currenttree;
naym *taxname=NULL, *speciesname=NULL, *psname=NULL;
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

  qvec = (int*)malloc(nquarts*sizeof(int));

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
      ppMatrix[i] = (double*)malloc(ntaxa*sizeof(double));
      if (ppMatrix[i]==NULL)
        {
          printf("Can't memalloc ppMatrix[%d]\n",i);
          return 0;
        }
    }

  ppBase_full = (int**)malloc(nseq*sizeof(int*));
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

  ppBase = (int**)malloc(nseq*sizeof(int*));
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
  for (i=0; i<3; i++)
   {
     ppNodeChildren[i] = (int*)malloc(ntaxa*sizeof(int));
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
     ppNodeChildrenLeftQuart[i] = (int*)malloc(ntaxa*sizeof(int));
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
     ppNodeChildrenRightQuart[i] = (int*)malloc(ntaxa*sizeof(int));
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
      fscanf(infile, "%s %s", psname[i],ppDNASeq[i]);
      if (verbose==1) printf("Read sequence %d with name %s\n",i,psname[i]);

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

  int i,j,k,s,m=0;
  int num_no_gaps=0;

  /* to facilitate bootstrapping later, put the constant A site pattern in the first position */
  /* in ppBase_unique; but don't add a count so that the total count will be correct          */
  /* for (i=0;i<nseq; i++) ppBase_unique[i][0] = 0; -- no -- no longer doing this             */

  for (j=0; j<nsite; j++) {

    if ((ppBase[0][j]!=99 && include_gaps==1) || (ppBase[0][j]!=99 && include_gaps==0 && ppBase[0][j]!=4)) {

      site_counter[m]=1;

      for (i=0; i<nseq; i++) ppBase_unique[i][m] = ppBase[i][j];

      ppBase[0][j]=99;

      for (k=j+1; k<nsite; k++) {

	i=0;
	
	while (i< nseq && ppBase_unique[i][m] == ppBase[i][k]) i++;

	if (i==nseq) {

	  site_counter[m]++;
	  ppBase[0][k]=99;

	}

	if (i<nseq && ppBase[i][k]==4 && include_gaps==0) ppBase[0][k]=99;

	if (i<nseq && include_gaps==0) {
	  
	  for (s=0; s<nseq; s++) {

	    if (ppBase[s][k]==4) ppBase[0][k]=99;

	  }

	}

      }

      m++;

    }

  }

  num_unique=m;
  printf("Finding unique site patterns ...");
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
  int nboot, ntquarts;
  int table_initialized = 0;
  naym tempname1, tempname2;
  float temp_rate;
  double totaltime;

  /* --- defaults --- */
  const char *settings_file   = "settings";
  const char *data_file       = "data.phy";
  const char *tree_file       = "treefile.tre";
  const char *out_file        = "outtree.tre";
  const char *picltrees_file  = "picltrees.tre";
  const char *results_file    = "results";
  const char *bootdata_file   = "boots.dat";


  /* --- override via argv --- */
  if (argc > 1) settings_file = argv[1];
  if (argc > 2) data_file     = argv[2];
  if (argc >3) tree_file      = argv[3];
  if (argc > 4) out_file      = argv[4];
  if (argc > 5) picltrees_file= argv[5];
  if (argc > 6) results_file  = argv[6];
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
  printf("\t\t\t\t Beta Version  \n");
  printf("\t\t\t         May 2026 \n");
  printf("----------------------------------------------------------------------------\n\n");

  // read from settings file
  set = fopen(settings_file,"r");
  
  /* Model */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&model);
  printf("Model: %d\n",model);
  
  /* Gaps */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&include_gaps);
  printf("Gaps: %d\n",include_gaps);

  /* Bootstrap */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&nboot);
  printf("Bootstrap: %d\n",nboot);

  /* Theta */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%f",&theta);
  printf("Theta: %f\n",theta);

  /* Lambda */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%f",&lambda);
  printf("Lambda: %f\n",lambda);
	
  /* Rate_param */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%f",&ratepar);
  printf("Rate_param: %f\n",ratepar);

  /* Random_tree */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&random_tree);
  printf("Random_tree: %d\n",random_tree);

  /* Opt_bl */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&anneal_bl);
  printf("Opt_bl: %d\n",anneal_bl);

  /* User_bl */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&user_bl);
  printf("User_bl: %d\n",user_bl);

  /* Num_opt */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&max_it_bl);
  printf("Num_Opt: %d\n",max_it_bl);

  /* Seed1 */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&seedj);
  printf("Seed1: %d\n",seedj);

  /* Seed2 */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&seedk);
  printf("Seed2: %d\n",seedk);

  /* Num_cat */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&ncat);
  printf("Num_cat: %d\n",ncat);

  /* Tree_search */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&anneal);
  printf("Tree_search: %d\n",anneal);

  /* Max_iter */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&max_it);
  printf("Max_iter: %d\n",max_it);
 
  /* Mult_iter */                                                                    
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&mult_iter);
  printf("Mult_iter: %d\n",mult_iter);
  //mult_iter = 3;

  /* prob_bound */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%lf",&prob_bound);
  printf("Prob_bound: %f\n",prob_bound);
  //prob_bound = 0.05;

  /* test_increment */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&test_increment);
  printf("Test_increment: %d\n",test_increment);
  //test_increment = 250; 

  /* b1opt */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%lf",&b1opt);
  printf("Target annealing slope: %f\n",b1opt);
  //b1opt = -0.01;
  
  /* Beta */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%f",&beta);
  printf("Beta: %f\n",beta);

  /* Verbose */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&verbose);
  printf("Verbose: %d\n\n",verbose);

  if (model == 5){ 
	  theta = NAN;
  	  
  }
  if (b1opt>0) printf("The target annealing slope (b1opt) must be negative; adaptive annealing is disabled.\n\n");

  fscanf(set,"%d",&ntaxa);
  setall((long)seedj,(long)seedk); /* Set seeds for random number generator */
  nquarts = binomial(ntaxa,4);
  printf("There are %d total species. This corresponds to %d species-level quartets.\n",ntaxa,nquarts);
  data = fopen(data_file,"r");
  fscanf(data,"%d %d\n",&nseq,&nsite);
  ntquarts = binomial(nseq,4);
  printf("There are %d total sequences. This corresponds to %d sequence-level quartets.\n\n",nseq,ntquarts);

  // allocate memory
  memcheck = MemAlloc();
  for (i=0; i<ntaxa+1; i++) seq_counter[i] = 0;

  taxname = (naym*)malloc(ntaxa*sizeof(naym));
  for (i=0; i<ntaxa; i++) fscanf(set,"%s",taxname[i]);

  // map sequences to species using info from settings file
  for (i=0; i<nseq; i++) {
        fscanf(set,"%s %s",tempname1,tempname2);
        j=0; while (strcmp(tempname1,taxname[j]) != 0) j++;
        ppSp_assign[j][seq_counter[j]] = i;
        seq_counter[j]++;
  }
  
  // checks that assignments were read in correctly
  for (j=0; j<ntaxa; j++) printf("Species %d with name %s has %d sequences\n",j,taxname[j],seq_counter[j]);
  printf("\n");
  fflush(0);

  printf("The species assignments are:\n");
  for (j=0; j<ntaxa; j++) {
	printf("Species %d: %s\n",j,taxname[j]);
	for (i=0; i<seq_counter[j]; i++) {
		printf("%d ",ppSp_assign[j][i]);
	}
	printf("\n");
  }
  printf("\n");

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
	if(model != 5){
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
  	printf("Tree succesfully read from file\n\n");        
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
  make_indmat();
  for (i=1; i<2*ntaxa; i++) {
	parents[i] = find_parent(i);
 	parents_temp[i] = parents[i];
  }

  /* test -- getting estimates of branch lengths via computation -- feature not yet implemented */
  //AllMyChildrenLeft(8);
  //AllMyChildrenRight(8);
  //AllMyChildren(8);
  //GetMomentEstimate(8);
  //for (i=0; i<3; i++) for (j=0;j<ntaxa;j++) printf("%d ",ppNodeChildren[i][j]);
  /* end test */

  /* Compute the composite likelihood for the current tree */
  if (model == 5) {
	  ComputeAandS_popvar(lambda);
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
		printf("Optimizing branch lengths using uphill search ...\n\n");
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
		printf("Optimizing branch lengths using simulating annealing ...\n\n");
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

	if (model == 1) printf("the composite likelihood of the tree is %f\n",GetCompLik());
  	else if (model == 2) printf("the composite likelihood of the tree is %f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);
  	else if (model == 3) printf("the composite likelihood of the tree is %f\n",GetCompLik_msnp());
	else if (model == 4) printf("the composite likelihood of the tree is %f\n",GetCompLik_genetree());
	else if (model == 5) printf("the composite likelihood of the tree is %f\n",GetCompLik_popvar());
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
  	}
  	for (i=0; i<2*ntaxa+1; i++) TimeVec[i] = TimeVec_best[i];
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
	printf("Results have been written to file outtree.tre -- exiting.\n\n");
	//exit(1);

  }
  else { printf("Tree_search must be 0, 1, or 2. Exiting.\n\n"); exit(1); }
  /* main analyses done */


  /* print branch lengths to file */
  if (verbose==1) {
  	res = fopen(results_file,"w");
	if (model == 5){
		for (i=1; i<ntaxa; i++) fprintf(res,"%f ",TimeVec[ntaxa+i]);
  	fclose(res);
	}
	  else{
		  for (i=1; i<ntaxa; i++) fprintf(res,"%f ",TimeVec[ntaxa+i]/theta);
  			fprintf(res,"%f \n",ratepar);
  			fclose(res);
	  }
  	
  }
  /* done print branch lengths to file */

  /* print only Newick string to file */
	pt = fopen(picltrees_file,"a");
	write_species_tree_out_file(ntaxa+1,ntaxa+1);
	fprintf(pt,";\n");
	fclose(pt);
  /* done print only Newick string to file */


  /* bootstrapping */
  if(nboot>0) {
	if (anneal_bl == 0) { printf("Please set Opt_bl = 1 to carry out bootstrapping -- exiting.\n\n"); exit(1); }
	else boot_times(nboot,bootdata_file);
  }

  printf("\n\nAnalysis complete - exiting.\n\n");

  /* useful functions -- not accessible to users through settings */
    //remove_CONSTANT();
    //remove_CONSTANT_MSNP();
    //print_PHYLIP();

}
