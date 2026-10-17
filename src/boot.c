/*************************************************************/
/*** Functions to carry out the bootstrap                  ***/
/*************************************************************/


void bootstrap(int nrep, double init_beta){

  int i, j, k, l, align_length;
  int nbatch, leftover, temptax;
  long *bs_sample, *subsample, *lsubsample;
  float *obs_site_pat_freqs;
  double lastbatch, **input_freqs, **last_input_freqs;
  float *finput_freqs, *flinput_freqs;


  printf("Beginning bootstrapping ....\n");

  /*******************/
  /* allocate memory */
  obs_site_pat_freqs = (float*)malloc((num_unique)*sizeof(float));
  if (obs_site_pat_freqs==NULL)
    {
      printf("     Can't memalloc obs_site_pat_freqs.\n");
    } 
  
  bs_sample = (long*)malloc((num_unique+1)*sizeof(long));
  if (bs_sample==NULL)
    {
      printf("     Can't memalloc bs_sample.\n");
    }

  subsample = (long*)malloc(200*sizeof(long));
  if (subsample==NULL)
   {
     printf("     Can't memalloc subsample.\n");
   } 

  finput_freqs = (float*)malloc(199*sizeof(float));
  if (finput_freqs==NULL)
   {
     printf("     Can't memalloc finput_freqs.\n");
   }

  /* done memalloc */
  /*****************/




  /*****************************************************/
  /* set length and initialize vector of probabilities */
  if (include_gaps==1) align_length = nsite;
  else {
         num_no_gaps = 0;
         for (i=0; i<num_unique; i++) num_no_gaps += site_counter[i];
         align_length = num_no_gaps;
  }

  /* When the number of unique sites is large, we need to break the multinomial categories up in */
  /* order to draw bootstrap samples - see Malefaki and Iliopoulos, CSDA 51: 5471-5476, 2007     */
  /* Do those in batchs of 200 after some testing with genmul                                    */
  /* The last batch will include 200 + whatever is leftover                                      */

  /* An annoying thing is that the genmul routine wants float precision for the input probs      */
  /* but when we normalize in each batch, we need double precision - so there's an extra step    */
  /* below to transfer to float                                                                  */

  nbatch = floor(num_unique/200);
  leftover = num_unique % 200;

  /* Memalloc for specific sequence features */
  if (nbatch>1) {
	input_freqs = (double**)malloc((nbatch-1)*sizeof(double*));
	for (k=0; k<nbatch-1; k++) input_freqs[k] = (double*)malloc(200*sizeof(double));
  }
  else {
	input_freqs = (double**)malloc(1*sizeof(double*));
	input_freqs[0] = (double*)malloc((200+leftover)*sizeof(double));
  }
  last_input_freqs = (double**)malloc(1*sizeof(double*));
  last_input_freqs[0] = (double*)malloc((200+leftover)*sizeof(double));

  double batch_sums[nbatch-1];
  float fbatch_sums[nbatch-1];
  long *batch_totals;

  if (nbatch>1) batch_totals = (long*)malloc((nbatch)*sizeof(long));
  else batch_totals = (long*)malloc(sizeof(long));
  if (batch_totals==NULL)
    {
      printf("     Can't memalloc batch_totals.\n");
    }
  
  flinput_freqs = (float*)malloc((200+leftover-1)*sizeof(float));
  if (flinput_freqs==NULL)
   {
     printf("     Can't memalloc flinput_freqs.\n");
   }

  lsubsample = (long*)malloc((200+leftover)*sizeof(long));
  if (lsubsample==NULL)
   {
     printf("     Can't memalloc lsubsample.\n");
   }
  /* Done memalloc */



  // get site pattern frequencies
  for (i=0; i<num_unique; i++) obs_site_pat_freqs[i] = (float)site_counter[i]/align_length; 

  // find probabilities in each batch of 200
  lastbatch = 1.0;
  if (nbatch > 0) {
  	for (k=0; k<nbatch-1; k++) {
		batch_sums[k]=0.0;
		for (l=0; l<200; l++) batch_sums[k] += obs_site_pat_freqs[200*k+l]; 
		lastbatch = lastbatch - batch_sums[k];
  	}
  }

  for (k=0; k<nbatch-1; k++) {
	fbatch_sums[k] = (float)batch_sums[k];
	for (l=0; l<200; l++) {
		input_freqs[k][l] = obs_site_pat_freqs[200*k+l]/batch_sums[k];
	}
  }

  // find probabilities for the last batch
  double tsum=0.0;
  for (l=0; l<200+leftover-1; l++) {
	last_input_freqs[0][l] = (double)obs_site_pat_freqs[200*(nbatch-1)+l]/(double)(lastbatch);
 	tsum += last_input_freqs[0][l];
	}
  last_input_freqs[0][200+leftover-1] = 1.0-tsum;

  /* done init */
  /*****************************************************/



  /*****************/
  /* now bootstrap */

  for (j=0; j<nrep; j++) {
  
        printf("Bootstrap rep: %d ",j+1);

	/* first generate number of sites in each batch */
	if (nbatch>0) genmul((long)align_length,fbatch_sums,nbatch,batch_totals);
	else batch_totals[0] = align_length;

	//printf("Batch totals: "); for (k=0; k<nbatch; k++) { printf("%ld ",batch_totals[k]);} printf("\n");

	long mycount;

	/* Then generate the specific sites in each batch */
	/* Use the first category as the missing one, to deal with float precision issues */
	/* Specifically, in genmul(), the sum of the category probabilities of the first  */
	/* n-1 categories must be less than 0.99999F. With most site patterns being       */
	/* and in general the site patterns being ordered in decreasing frequency, this   */
	/* is commonly violated. Putting the first site pattern "at the end" takes care   */
	/* of this issue.								  */


	for (k=0; k<nbatch-1; k++) {

		for (l=1; l<200; l++) finput_freqs[l-1] = (float)input_freqs[k][l];
        	genmul(batch_totals[k],finput_freqs,200,subsample);

		for (l=1; l<200; l++) site_counter[200*k+l] = subsample[l-1];
		site_counter[200*k] = subsample[199];

	}

	/* The last batch is done separately bc it has a different length */
	/* Since the very last pattern has the smallest probabaility,     */
        /* which can cause problems in genmul, use the first category     */
        /* as the missing one                                             */

	//printf("need to sample %ld more\n",batch_totals[nbatch-1]);

        for (l=1; l<200+leftover; l++) { 
		flinput_freqs[l-1] = (float)last_input_freqs[0][l];
	}

	genmul(batch_totals[nbatch-1],flinput_freqs,200+leftover,lsubsample);
	for (l=1; l<200+leftover; l++) site_counter[200*(nbatch-1)+l] = lsubsample[l-1];
	site_counter[200*(nbatch-1)] = lsubsample[200+leftover-1];

	int bootcount=0;
	for (i=0; i<num_unique; i++) bootcount+=site_counter[i];

	/* The site_counter vector is resampled for each bootstrap   */
        /* replicate, so the cached per-species quartet site pattern */
        /* counts in complik.c must be deleted and re-computed       */

	InvalidateQuartetSitePatternCache();

	/* The bootstrap sample is generated, now do estimation */

	// copy starting tree to besttree before beginning annealing
	// copy the tree estimated from the initial annealing to ppTwoRow and TimeVec
        for (i=0; i<ntaxa; i++) {
		ppTwoRow[0][i] = ppTwoRow_best_orig[0][i];
		ppTwoRow[1][i] = ppTwoRow_best_orig[1][i];
		ppTwoRow_temp[0][i] = ppTwoRow_best_orig[0][i];
		ppTwoRow_temp[1][i] = ppTwoRow_best_orig[1][i];
                ppTwoRow_best[0][i] = ppTwoRow[0][i];
                ppTwoRow_best[1][i] = ppTwoRow[1][i];
        }
        for (i=0; i<2*ntaxa+1; i++) {
		TimeVec[i] = TimeVec_best_orig[i];
		TimeVec_temp[i] = TimeVec_best_orig[i];
		TimeVec_best[i] = TimeVec[i];
	}
	make_indmat();
  	for (i=1; i<2*ntaxa; i++) {
        	parents[i] = find_parent_scan(i);
        	parents_temp[i] = parents[i];
  	}

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
  	//printf("The composite likelihood of the initial tree is %f\n\n",curr_anneal_lik);

	/* restore beta to starting value and anneal */

	beta = init_beta;
        if (model == 1) anneal_full();
	else if (model == 2) anneal_ratevar();
        else if (model == 3) anneal_msnp();
	else if (model == 4) anneal_genetree();
	else if (model == 5) anneal_full_popvar();

	/* Re-order besttree and write to file */
	for (i=0; i<ntaxa; i++) {
                if (ppTwoRow_best[0][i]>ppTwoRow_best[1][i]) {
                temptax = ppTwoRow_best[0][i];
                ppTwoRow_best[0][i] = ppTwoRow_best[1][i];
                ppTwoRow_best[1][i] = temptax;
                }
        }
 	write_species_tree_boot(ntaxa+1,ntaxa+1);
        fprintf(boot,";\n");
	fflush(0);

   }
   /* done bootstrap */
   /******************/

   //fclose(boot);
   printf("\n\n Done bootstrapping. Results have been written to boots.dat.\n\n");

}
