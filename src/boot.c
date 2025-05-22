/*************************************************************/
/*** Functions to carry out the bootstrap                  ***/
/*************************************************************/


void boot_times(int nrep){

  int i, j, k, l, align_length;
  int nbatch, leftover;
  long *bs_sample, *subsample, *lsubsample;
  float *obs_site_pat_freqs;
  double lastbatch, **input_freqs;
  float *finput_freqs, *flinput_freqs;
  FILE *boot;
	
  boot = fopen("boot.dat","w");
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

  //printf("Alignment length is %d and num_unique is %d and num_no_gaps is %d\n",align_length,num_unique,num_no_gaps);

  /* When the number of unique sites is large, we need to break the multinomial categories up in */
  /* order to draw bootstrap samples - see Malefaki and Iliopoulos, CSDA 51: 5471-5476, 2007     */
  /* Do those in batchs of 200 after some testing with genmul                                    */

  /* An annoying thing is that the genmul routine wants float precision for the input probs      */
  /* but when we normalize in each batch, we need double precision - so there's an extra step    */
  /* below to transfer to float                                                                  */

  nbatch = floor(num_unique/200);
  leftover = num_unique % 200;

  /* Memalloc for specific sequence features */
  input_freqs = (double**)malloc((nbatch+1)*sizeof(double*));
  for (k=0; k<nbatch+1; k++) input_freqs[k] = (double*)malloc(200*sizeof(double));

  double batch_sums[nbatch];
  float fbatch_sums[nbatch];
  //long *batch_totals[nbatch+1];
  long *batch_totals;

  batch_totals = (long*)malloc((nbatch+1)*sizeof(long));
  if (batch_totals==NULL)
    {
      printf("     Can't memalloc batch_totals.\n");
    }
  
  flinput_freqs = (float*)malloc((leftover-1)*sizeof(float));
  if (flinput_freqs==NULL)
   {
     printf("     Can't memalloc flinput_freqs.\n");
   }

  lsubsample = (long*)malloc(leftover*sizeof(long));
  if (lsubsample==NULL)
   {
     printf("     Can't memalloc lsubsample.\n");
   }
  /* Done memalloc */

  //printf("The number of batches is %d and the leftover is %d\n",nbatch,leftover);
 
  // get site pattern frequencies
  for (i=0; i<num_unique; i++) obs_site_pat_freqs[i] = (float)site_counter[i]/align_length; 

  // temp check
  //  printf("The observed pattern frequencies are:\n");
  // double sum = 0.0;
  // for (i=0; i<num_unique; i++) sum += obs_site_pat_freqs[i];
  // printf("sum is %f\n",sum);	 
  //for (i=0; i<num_unique; i++) printf("%f ",obs_site_pat_freqs[i]);

  // find probability in each batch of 200
  lastbatch = 1.0;
  if (nbatch > 0) {
  	for (k=0; k<nbatch; k++) {
		batch_sums[k]=0.0;
		for (l=0; l<200; l++) batch_sums[k] += obs_site_pat_freqs[200*k+l]; 
		lastbatch = lastbatch - batch_sums[k];
  	}
  }

  double check_probs;
  for (k=0; k<nbatch; k++) {
	fbatch_sums[k] = (float)batch_sums[k];
	check_probs=0.0;
	for (l=0; l<200; l++) {
		input_freqs[k][l] = obs_site_pat_freqs[200*k+l]/batch_sums[k];
		check_probs += input_freqs[k][l];
	}
	//printf("The total prob for batch %d is %f\n",k,check_probs);
  }


  double tsum=0.0;
  for (l=0; l<leftover-1; l++) {
	input_freqs[nbatch][l] = (double)obs_site_pat_freqs[200*nbatch+l]/(double)lastbatch;
 	tsum += input_freqs[nbatch][l];
	}
  input_freqs[nbatch][leftover-1] = 1.0-tsum;
  for (l=leftover; l<200; l++) input_freqs[nbatch][l] = 0.0;

  double temp_check=0.0;
  double minn= 1.0;
  for (l=0; l<leftover-1; l++) {
	temp_check+=input_freqs[nbatch][l];
	if (input_freqs[nbatch][l]<minn) minn = input_freqs[nbatch][l];
	}
  //printf("temp_check %f,total is %f\n",temp_check,temp_check+input_freqs[nbatch][leftover-1]);
  //printf("minn is %f %f\n",minn,(float)minn);

  /* done init */
  /*****************************************************/



  /*****************/
  /* now bootstrap */

  printf("\n Bootstrap rep: ");
  for (j=0; j<nrep; j++) {
  
        printf("%d ",j);

	//for (k=0; k<nbatch; k++) printf("%f ",batch_sums[k]);
	//printf("Align length is %ld\n",(long)align_length);

	/* first generate number of sites in each batch */
	if (nbatch>0) genmul((long)align_length,fbatch_sums,nbatch+1,batch_totals);
	else batch_totals[0] = align_length;

	long mycount;

	/* Then generate the specific sites in each batch */
	for (k=0; k<nbatch; k++) {

		//printf("%ld ",batch_totals[k]);

		for (l=0; l<199; l++) finput_freqs[l] = (float)input_freqs[k][l];

        	genmul(batch_totals[k],finput_freqs,200,subsample);

		mycount=0;
		for (l=0; l<200; l++) mycount += subsample[l];

  		//printf("For batch %d, total number sampled is %ld\n",k,mycount);

		for (l=0; l<200; l++) site_counter[200*k+l] = subsample[l];

	}

	/* The last batch is done separately bc it has a different length */
	/* Since the very last pattern has the smallest probabaility,     */
        /* which can cause problems in genmul, use the first category     */
        /* as the missing one                                             */

	//printf("need to sample %ld more; nbatch is %d\n",batch_totals[nbatch],nbatch);

	float check_sum = 0.0;
	double my_check_sum = 0.0;
        for (l=1; l<leftover; l++) { 
		flinput_freqs[l-1] = (float)input_freqs[nbatch][l];
		check_sum += flinput_freqs[l-1];
		my_check_sum += input_freqs[nbatch][l];
	}
        //printf("checksum %f, mychecksum is %f\n",check_sum,my_check_sum);
	genmul(batch_totals[nbatch],flinput_freqs,leftover,lsubsample);
	for (l=1; l<leftover; l++) site_counter[200*nbatch+l] = lsubsample[l-1];
	site_counter[200*nbatch] = lsubsample[leftover-1];

	int bootcount=0;
	for (i=0; i<num_unique; i++) bootcount+=site_counter[i];
	//printf("Total number of sampled sites is %d\n",bootcount);

	/* bootstrap sample is generated, now do estimation */
	for (i=1; i<2*ntaxa+1; i++) TimeVec[i] = TimeVec_init[i];

        if (anneal == 0) {
        	if (model == 1) bl_uphill_full();
		else if (model == 2) bl_uphill_ratevar();
        	else if (model == 3) bl_uphill_msnp();
        }
  	else {
        	if (model == 1) bl_anneal_full();  
		else if (model == 2) bl_anneal_ratevar();
        	else if (model == 3) bl_anneal_msnp();
        }

        for (i=0; i<ntaxa-1; i++) fprintf(boot,"%1.12f ",TimeVec[i+ntaxa+1]/theta);
        fprintf(boot,"\n");
	fflush(0);

   }
   /* done bootstrap */
   /******************/

   fclose(boot);
   printf("\n\n Done bootstrapping. Results have been written to boot.dat.\n\n");

}
