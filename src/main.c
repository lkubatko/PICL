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
#include "nodeopt.c"
#include "anneal.c"
#include "boot.c"
#include "opt.c"
#include "trbldg.c"
#include "trbldg_ratevar.c"
#include "trbldg_msnp.c"

// ntaxa = number of species
// nseq = number of individuals
int ntaxa, nspecies, nseq, nsite, nquarts, num_unique_quarts, num_unique, include_gaps, num_no_gaps, verbose, model, ncat, random_tree;
int anneal, anneal_bl, user_bl, max_it, max_it_bl, seedj, seedk;
int *parents, *parents_temp, *ppTwoRow[2], *ppTwoRow_temp[2], *ppTwoRow_best[2], *ppTwoRowQuart[2], *filled_ind, *seq_counter, *qvec;
int **ppBase_full, **ppBase, **ppBase_unique, *site_counter, **ppSp_assign, **ppNodeChildren, **ppNodeChildrenLeftQuart, **ppNodeChildrenRightQuart;
double ci, max_cl;
float theta, beta, mu, ratepar, invpar;
double *TimeVec, *TimeVec_temp, *TimeVec_init, *TimeVec_best, *TimeVecQuart, *rvals, **ppLengthMat, **ppMatrix;
double smat[10][10],amat[12][12];
FILE *out;

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

  site_counter = (int*)malloc((nsite+1)*sizeof(int));
  if (site_counter==NULL)
    {
      printf("     Can't memalloc site_counter.\n");
      return 0;
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

  for (i=0;i<nseq; i++) ppBase_unique[i][0] = 0;

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
  if (include_gaps == 0) printf("Sites with a gap in at least one taxa have been excluded - the total number of sites used in the computations is %d\n\n",num_no_gaps);

}


int main() {

  FILE *set, *tf, *data, *res;
  char keyword;
  int memcheck, i, j, k, l;
  int i2, j2, k2, l2;
  int nboot;
  naym tempname1, tempname2;
  float temp_rate;
  double totaltime;

  // print opening message
  printf("\n\n----------------------------------------------------------------------------\n");
  printf("\t PICL: Phylogenetic Inference with Composite Likelihood \n");
  printf("\t\t\t Very Beta Version  \n");
  printf("\t\t\t      May 2025 \n");
  printf("----------------------------------------------------------------------------\n\n");

  // read from settings file
  set = fopen("settings","r");
  
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

  /* Num_iter */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&max_it);
  printf("Num_iter: %d\n",max_it);

  /* Beta */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%f",&beta);
  printf("Beta: %f\n",beta);

  /* Verbose */
  keyword = 0; while (keyword != 58) keyword = fgetc(set); fgetc(set); fscanf(set,"%d",&verbose);
  printf("Verbose: %d\n\n",verbose);


  fscanf(set,"%d",&ntaxa);
  setall((long)seedj,(long)seedk); /* Set seeds for random number generator */
  nquarts = binomial(ntaxa,4);
  printf("There are %d total sequences. This corresponds to %d species-level quartets.\n\n",ntaxa,nquarts);
  data = fopen("data.phy","r");
  fscanf(data,"%d %d\n",&nseq,&nsite);

  // allocate memory
  memcheck = MemAlloc();

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
  else if (model==2) printf("Estimation will be carried out using the CIS/multlilocus (MSC-JC69 + G model)\n");
  else if (model==3) printf("Estimation will be carried out using the SNP model\n");
  else { printf("Model must be 1, 2, or 3. Exiting.\n\n"); exit(1); }

  /******* TO DO: need to free extra memory from storing taxon map *******/

  // read in data
  getSequenceData(data);
  fclose(data);
  printf("\nSequence alignment with %d lineages and %d sites has been read successfully\n",nseq,nsite);
  unique_sites();

  // set a tree
  if (random_tree==0) { // read in user species tree
  	for (i=0; i<2*ntaxa; i++) TimeVec[i] = 0.0;
  	tf = fopen("treefile.tre","r");
	if (tf==NULL) { printf("File treefile.tre not found. Exiting.\n\n"); exit(1);}
  	temp_rate = ReadTree(tf);
  	totaltime = FindTotalTime();
  	CalcTimeVec(totaltime,1.0);
	
  	// Transfer from temp and re-scale with theta if needed

  	for (i=0; i<ntaxa-1; i++) {
		ppTwoRow[0][i] = ppTwoRow_temp[0][i];
	 	ppTwoRow[1][i] = ppTwoRow_temp[1][i];
  	}
  	for (i=0; i<2*ntaxa+1; i++) {
		TimeVec[i] = theta*TimeVec_temp[i];
		TimeVec_init[i] = theta*TimeVec_temp[i];
		TimeVec_temp[i] = TimeVec[i];
  	}

  	if (verbose == 1 ) {
    		printf("The tree read from the input file is (coalescent units, mutation units):\n");
    		for (i=0; i<ntaxa-1; i++) {
  		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
  		printf("%f %f\n ",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
  		}
		printf("\n");
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
                printf("The random tree generated is (coalescent units, mutation units):\n");
                for (i=0; i<ntaxa-1; i++) {
                printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
                printf("%f %f\n ",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
        	}
	printf("\n");
	}
   }
   else {  // print an error and exit
	printf("Random_tree must be either 0 or 1. Exiting.\n\n");
	exit(1);
   }

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
  ComputeAandS(theta);
  for (i=0; i<nquarts+1; i++) qvec[i]=0; 
  if (model == 1) printf("The composite likelihood of the initial tree is %f\n",GetCompLik());
  else if (model == 2) printf("The composite likelihood of the initial tree is %f\n",GetCompLik_ratevar());
  else if (model == 3) printf("The composite likelihood of the initial tree is %f\n",GetCompLik_msnp());
  printf("\n");

  /*************  Set-up complete ***************/


  /*************  Run the appropriate analysis **************/

  // Evalute current tree only
  // if anneal = 0
  if (anneal==0) {
	if (anneal_bl==0) { printf("Analysis complete -- exiting.\n\n"); exit(1); }
	else if (anneal_bl==1) { // uphill bl search
		if (user_bl == 0) assign_times(); //use equally-spaced starting times
		if (model == 1) bl_uphill_full();
        	else if (model == 2) {
                	bl_uphill_ratevar();                                                                                               
                	OptAlpha();
                	bl_uphill_ratevar();
        	}
        	else if (model == 3) bl_uphill_msnp();

	}
	else if (anneal_bl==2) { // simulated annealing bl search
		if (user_bl == 0) assign_times(); //use equally-spaced starting times
		if (model == 1) bl_anneal_full();
        	else if (model == 2) {
                	bl_anneal_ratevar();
                	OptAlpha();
                	bl_anneal_ratevar();
        	}
        	else if (model == 3) bl_anneal_msnp();
	}
	else if (anneal_bl==3) { // numerical optimization of branch lengths
		printf("Numerical optimization of branch lengths is not yet implemented. Exiting.\n\n"); exit(1);
	}
	else { printf("Opt_bl must be 0, 1, 2, or 3. Exiting.\n\n"); exit(1); }

        if (verbose == 1 ) {
    		printf("The tree after branch length optimization is (coalescent units, mutation units):\n");
    		for (i=0; i<ntaxa-1; i++) {
        		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        		printf("%f %f\n",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
		}
	}
	printf("After branch length optimization, ");
	if (model == 1) printf("the composite likelihood of the tree is %f\n",GetCompLik());
  	else if (model == 2) printf("the composite likelihood of the tree is %f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);
  	else if (model == 3) printf("the composite likelihood of the tree is %f\n",GetCompLik_msnp());
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

        out = fopen("outtree.tre","w");
        fprintf(out,"Current tree in mutation units: \n");
        write_species_tree_out(ntaxa+1,ntaxa+1);
	fprintf(out,";\n\n");
	fprintf(out,"Current tree in coalescent units: \n");
        for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta;
        write_species_tree_out(ntaxa+1,ntaxa+1);
        fprintf(out,";\n\n");
        for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]*theta;
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
	if (model == 1) {
		bl_uphill_full();
		anneal_full();
		if (anneal_bl==1) { printf("optimizing bl%f\n",GetCompLik()); bl_uphill_full(); printf("fdone %f\n\n",GetCompLik()); }
		else if (anneal_bl==2) bl_anneal_full();
	} 
	//else if (model == 2) anneal_ratevar();
	//else if (model == 3) anneal_msnp();

	printf(" done\n\n");
	if (verbose==1) {
		printf("After search, tree is \n");
  		for (i=0; i<ntaxa-1; i++) {
        		printf("%d %d ",ppTwoRow[0][i],ppTwoRow[1][i]);
        		printf("%f %f\n",TimeVec[i+ntaxa+1]/theta,TimeVec[i+ntaxa+1]);
  		}
		if (model == 1) printf("The composite likelihood of this tree is %f\n",GetCompLik());
                else if (model == 2) printf("The composite likelihood of this tree is %f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);
                else if (model == 3) printf("The composite likelihood of this tree is %f\n",GetCompLik_msnp());
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

       	out = fopen("outtree.tre","w");
  	fprintf(out,"Current tree in mutation units: \n");
  	write_species_tree_out(ntaxa+1,ntaxa+1);
	fprintf(out,";\n\n");
	fprintf(out,"Current tree in coalescent units: \n");
	for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta;
	write_species_tree_out(ntaxa+1,ntaxa+1);
  	fprintf(out,";\n\n");   
	for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]*theta; 
  
  	// copy best tree to current tree and write to file outtree.tre
  	for (i=0; i<ntaxa-1; i++) {
         	ppTwoRow[0][i] = ppTwoRow_best[0][i];
         	ppTwoRow[1][i] = ppTwoRow_best[1][i];
  	}
  	for (i=0; i<2*ntaxa+1; i++) TimeVec[i] = TimeVec_best[i];

	printf("The composite likelihood of the best tree found by the algorithm is ");
	if (model == 1) printf("%f\n",GetCompLik());
        else if (model == 2) printf("%f with rate variation parameter %f\n",GetCompLik_ratevar(),ratepar);    
        else if (model == 3) printf("%f\n",GetCompLik_msnp());
        printf("\n\n");
         
   	fprintf(out,"Best tree found by the algorithm in mutation units: \n");
  	write_species_tree_out(ntaxa+1,ntaxa+1);
  	fprintf(out,";\n\n");
	fprintf(out,"Best tree found by the algorithm in coalescent units: \n");
        for (i=1;i<ntaxa; i++) TimeVec[ntaxa+i] = TimeVec[ntaxa+i]/theta; 
        write_species_tree_out(ntaxa+1,ntaxa+1);     
	fprintf(out,";\n\n");

	fclose(out);
	printf("Results have been written to file outtree.tre -- exiting.\n\n");
	exit(1);

  }
  else { printf("Tree_search must be 0, 1, or 2. Exiting.\n\n"); exit(1); }
  /* main analyses done */


  /* print branch lengths to file */
  if (verbose==1) {
  	res = fopen("results","w");
  	for (i=1; i<ntaxa; i++) fprintf(res,"%f ",TimeVec[ntaxa+i]/theta);
  	fprintf(res,"%f \n",ratepar);
  	fclose(res);
  }
  /* done print branch lengths to file */


  /* bootstrapping */
  if(nboot>0) {
	if (anneal_bl == 0) { printf("Please set Opt_bl = 1 to carry out bootstrapping -- exiting.\n\n"); exit(1); }
	else boot_times(nboot);
  }

  printf("\n\nAnalysis complete - exiting.\n\n");

  /* useful functions -- not accessible through settings */
  //  remove_CONSTANT();
  //  print_PHYLIP();

}
